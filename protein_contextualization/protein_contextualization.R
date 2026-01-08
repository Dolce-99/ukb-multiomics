library(prodente)
library(dplyr)
library(data.table)
library(ggplot2)

grab.sig = function(trait, covar, bonferroni){
    prote.feat = fread(paste0('mp_',trait,'_',covar,'_weighted_abs.txt')) %>% 
                  select(-V1) %>% filter(Block == "prote") # only select protein
    prote.feat = prote.feat %>% mutate(Imp = Imp_Raw*sign(Correlation), norm = scale(Imp, center=FALSE), protein_name = gsub("_INVN","", Feature)) %>% 
                    arrange(desc(abs(norm)))
    top.prote.feat = prote.feat %>% filter(abs(norm) > 1.96) %>% mutate(protein_name = tolower(protein_name)) %>% filter(protein_name != "batch")
    print(nrow(top.prote.feat))
    check_protein_overlap(top.prote.feat$protein_name, return_missing = T)
    results <- enrich_protein_characteristics(
                protein_foreground = top.prote.feat$protein_name, 
                factor_minimum_explained_variance = 0, 
                n_cores = 8)
    results = results %>% filter(population == "All")
    if (bonferroni == TRUE){
      results[, p_adjust := stats::p.adjust(pval, method = "bonferroni")]
      sig = results %>% filter(p_adjust < 0.05 & or > 1) %>% arrange(desc(or), p_adjust)
    } else {
      sig = results %>% filter(pval < 0.05 & or > 1) %>% arrange(desc(or), pval)
    }
    sig$output.label = paste0(sig$label, " (OR=", round(sig$or, 2), ", adjusted p-value=", sprintf("%.2e", sig$p_adjust), ")")
    return(sig)
}

res_bon = data.frame()
for (trait in c('MACE','Heart_Failure', 'CHD', 'PAD', 'Atrial_Fibrillation', 'T2_Diabetes',
                'Prostate_Cancer', 'Skin_Cancer', 'Breast_Cancer',
                'Dementia', 
                'Cataracts',  'Glaucoma',
                'COPD', 'Asthma', 
                'Renal_Disease', 'Liver_Disease', 'Fractures')){
      for (covarset in c("age_sex","ASCVD","PANEL")){
          sig_bon = grab.sig(trait, covar = covarset, bonferroni = TRUE)
          for (type in c("Drugs", "Socioeconomic", "Bone", "Diet", "Health", "Cardiovascular", "Blood cell counts", "pollution", "Basic demographics", "Ancestry", "Pulmonary", "Genetic", "Biomarker","Technical")){
              sig_bon_type = sig_bon %>% filter(category == !!type)
              res.dt_bon = data.table(Trait = trait, Covar = covarset, Type = type, label = paste(sig_bon_type$output.label,collapse = ", "))
              res.dt_bon[, label := sapply(label, toString)]
              res_bon = rbind(res_bon, res.dt_bon)
      }
  }
}


for (type in c("Drugs", "Socioeconomic", "Bone", "Diet", "Health", "Cardiovascular", "Blood cell counts", "pollution", "Basic demographics", "Ancestry", "Pulmonary", "Genetic", "Biomarker","Technical")){
    res_bon_type = res_bon %>% filter(Type == !!type) %>% select(-Type)
    write.table(as.data.frame(res_bon_type), paste0(type, '_enrichment_bon.txt'), sep = "\t", quote = F, row.names = F, col.names = T)
}
