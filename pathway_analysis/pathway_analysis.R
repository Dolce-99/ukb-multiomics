library(tidyverse)
library(clusterProfiler)
library(GO.db)
library(org.Hs.eg.db)
library(data.table)

backgroundgenes <- fread("backgroundgenes.csv")
bkgd.genes <- backgroundgenes[,2]
bkgd.entrez = bitr(bkgd.genes$V2, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
genes_bg = paste0("hsa:", bkgd.entrez$ENTREZID)

# process metabolites
metab_list <- list()
for (trait in c('MACE','Heart_Failure', 'CHD', 'PAD', 'Atrial_Fibrillation', 'T2_Diabetes',
                'Prostate_Cancer', 'Skin_Cancer', 'Breast_Cancer',
                'Dementia', 
                'Cataracts',  'Glaucoma',
                'COPD', 'Asthma', 
                'Renal_Disease', 'Liver_Disease', 'Fractures')){
    feat.metab = fread(paste0("/proj/yunligrp/users/djw/UKB_multiomics/revision1/02.new_features/mp_",trait,"_PANEL_weighted_abs.txt")) %>% 
                dplyr::select(-V1) %>% filter(Block == "meta") %>% mutate(Imp = Imp_Raw*sign(Correlation)) %>% 
                mutate(norm = scale(Imp), metab_name = Feature)
    top.metab = feat.metab %>% filter(abs(norm) > 1.96) %>% pull(metab_name)                
    metab_list[[trait]] <- top.metab
}
used.metab = unique(unlist(metab_list))

# manual mapping
proxy_map <- data.frame(
  Feature = c("Glucose", "Lactate", 
              "LA", "Omega.6", "PUFA", 
              "Omega.3", 
              "Total.C", "Clinical.LDL.C", "non.HDL.C", 
              "Total.TG", "Total.L", "VLDL.L", 
              "Total.FA", "SFA", "MUFA",
              "Albumin"),
  
  KEGG_ID = c("cpd:C00031", "cpd:C00186", 
              "cpd:C01595", "cpd:C01595", "cpd:C01595",
              "cpd:C06429",
              "cpd:C00187", "cpd:C00187", "cpd:C00187",
              "cpd:C00165", "cpd:C00165", "cpd:C00165",
              "cpd:C00249", "cpd:C00249", "cpd:C00249",
              "hsa:213")
)
metab_bg <- unique(proxy_map$KEGG_ID)
bkgd_universe = c(genes_bg, metab_bg)

link_gene <- keggLink("pathway", "hsa")
df_gene <- data.frame(
  term = sub("path:", "", unname(link_gene)), 
  gene = sub("path:", "", names(link_gene))  
)

link_cpd <- keggLink("pathway", "compound")
df_cpd <- data.frame(
  term = sub("path:map", "hsa", unname(link_cpd)), 
  gene = sub("cpd:", "cpd:", names(link_cpd)) 
)

term2gene_hybrid <- rbind(df_cpd, df_gene)

path_names <- keggList("pathway", "hsa")
term2name_hybrid <- data.frame(
  term = sub("path:", "", names(path_names)), 
  name = path_names,                     
  stringsAsFactors = FALSE
)

run_share_pathway = function(trait, map_df, universe_list, term2gene_db, term2name_db){
    optimized_universe <- intersect(universe_list, unique(term2gene_db$gene))
    feat.prote = fread(paste0("mp_",trait,"_PANEL_weighted_abs.txt")) %>% 
                dplyr::select(-V1) %>% filter(Block == "prote") %>% mutate(Imp = Imp_Raw*sign(Correlation)) %>% 
                mutate(norm = scale(Imp,center=FALSE), protein_name = gsub("_INVN","", Feature))
    top.prote = feat.prote %>% dplyr::filter(abs(norm) > 1.96) %>% pull(protein_name)
    genes = bitr(top.prote, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    prote.kegg = enricher(gene = paste0("hsa:", genes$ENTREZID),
                        universe = optimized_universe,
                        TERM2GENE = term2gene_db,
                        TERM2NAME = term2name_db,
                        pAdjustMethod = "fdr",
                        pvalueCutoff = 1)
    n_prote_input = length(top.prote)
    n_prote_mapped = length(intersect(paste0("hsa:", genes$ENTREZID), universe_list))

    feat.metab = fread(paste0("mp_",trait,"_PANEL_weighted_abs.txt")) %>% 
                dplyr::select(-V1) %>% filter(Block == "meta") %>% mutate(Imp = Imp_Raw*sign(Correlation)) %>% 
                mutate(norm = scale(Imp, center=FALSE), metab_name = Feature)
    top.metab = feat.metab %>% filter(abs(norm) > 1.96) %>% pull(metab_name)
    mapped_ids <- map_df %>% filter(Feature %in% top.metab) %>% pull(KEGG_ID) %>% unique()
    metab.kegg = enricher(gene = mapped_ids,
                              universe = optimized_universe,
                              TERM2GENE = term2gene_db,
                              TERM2NAME = term2name_db,
                              pvalueCutoff = 1,
                              pAdjustMethod = "fdr")
    n_metab_input = length(top.metab)
    n_metab_mapped = length(intersect(mapped_ids, universe_list))
    return(list(
        prote_res = prote.kegg, 
        metab_res = metab.kegg,
        counts = list(
            prote_in = n_prote_input,
            prote_used = n_prote_mapped,
            metab_in = n_metab_input,
            metab_used = n_metab_mapped
        )
    ))
}


summary_rows = list()
protein_results = list()
metab_results = list()
for (trait in c('MACE','Heart_Failure', 'CHD', 'PAD', 'Atrial_Fibrillation', 'T2_Diabetes',
                'Prostate_Cancer', 'Skin_Cancer', 'Breast_Cancer',
                'Dementia', 
                'Cataracts',  'Glaucoma',
                'COPD', 'Asthma', 
                'Renal_Disease', 'Liver_Disease', 'Fractures')){
    kegg.res = run_share_pathway(trait, proxy_map, bkgd_universe, term2gene_hybrid, term2name_hybrid)
    res_prote <- kegg.res$prote_res
    res_metab <- kegg.res$metab_res
    counts <- kegg.res$counts
    if (!is.null(res_prote) && nrow(as.data.frame(res_prote)) > 0) {
        df_prote <- as.data.frame(res_prote) %>% filter(pvalue < 0.05)
        n_prote <- nrow(df_prote)
        desc_prote <- df_prote$Description
    } else {
        n_prote <- 0
        desc_prote <- character(0)
    }
  
  # Metabolite processing
    if (!is.null(res_metab) && nrow(as.data.frame(res_metab)) > 0) {
        df_metab <- as.data.frame(res_metab) %>% filter(pvalue < 0.05)
        n_metab <- nrow(df_metab)
        desc_metab <- df_metab$Description
    } else {
        n_metab <- 0
        desc_metab <- character(0)
    }
  
    shared_terms <- intersect(desc_prote, desc_metab)
  
    clean_shared <- gsub(" - Homo sapiens \\(human\\)", "", shared_terms)
  
    if (length(clean_shared) > 0) {
        shared_string <- paste(clean_shared, collapse = ", ")
    } else {
        shared_string <- "" # Empty if no overlap
    }
  
     summary_rows[[trait]] <- data.frame(
        Trait = trait,
        Proteins_Inputted = counts$prote_in,
        Proteins_Used = counts$prote_used,
        Metabolites_Inputted = counts$metab_in,
        Metabolites_Used = counts$metab_used,
        
        # Result Columns
        N_Significant_Pathway_Protein = n_prote,
        N_Significant_Pathway_Metab = n_metab,
        Shared_Pathways = shared_string,
        stringsAsFactors = FALSE
    )

    df_prote = df_prote %>% mutate(Description = gsub(" - Homo sapiens \\(human\\)", "", Description)) %>% 
                mutate(Trait = trait) %>% dplyr::select(Trait, Description, GeneRatio, BgRatio, RichFactor, FoldEnrichment, zScore, pvalue, p.adjust, qvalue, geneID, Count)
    df_metab = df_metab %>% mutate(Description = gsub(" - Homo sapiens \\(human\\)", "", Description)) %>% 
                mutate(Trait = trait) %>% dplyr::select(Trait, Description, GeneRatio, BgRatio, RichFactor, FoldEnrichment, zScore, pvalue, p.adjust, qvalue, geneID, Count)
    protein_results[[trait]] = df_prote
    metab_results[[trait]] = df_metab
}

final_df <- do.call(rbind, summary_rows)
rownames(final_df) <- NULL
print(final_df)

protein_final_df = do.call(rbind, protein_results)
metab_final_df = do.call(rbind, metab_results)
rownames(protein_final_df) = NULL
rownames(metab_final_df) = NULL
write.table(protein_final_df, "protein_pathway_results.txt", sep="\t", quote=F, row.names = F)
write.table(metab_final_df, "metab_pathway_results.txt", sep="\t", quote=F, row.names=F)


