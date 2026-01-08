library(survival)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(survcomp)

covar.imp = fread("covariates_imp_new.tab")
tte = fread("ICD10_time_to_event.tab")
comb.df = tte %>% right_join(covar.imp,by="f.eid")
comb.df$baseline_visit = as.Date(comb.df$baseline_visit)
comb.df$DOB = as.Date(comb.df$DOB)
randomseeds <- c("rs42", "rs123", "rs456", "rs1024", "rs2048")
res.type = 'fix'

adj = function(covar.df, trait, covar.set.type){
    trait.date = paste0(trait,"_date")
    covar.df[[trait.date]] = as.Date(covar.df[[trait.date]])
    # for each trait, remove individuals with cases before visit
    covar.df = covar.df %>% filter(baseline_visit < !!sym(trait.date)) %>% mutate(event.age = as.numeric(!!sym(trait.date) - baseline_visit))
    
    # identify sex specific traits
    remove_sex = trait %in% c("Breast_Cancer","Prostate_Cancer")
    covar.set = switch(covar.set.type,
                age_sex = {
                    if(remove_sex) "Age + as.factor(center)" else "Age + as.factor(sex) + as.factor(center)"
                },
                ASCVD = {
                    if(remove_sex) 
                        "Age + as.factor(Current_Smoker) + as.factor(T2D) + Systolic_Blood_Pressure + Total_Cholesterol + HDL_Cholesterol + as.factor(BP_medication) + as.factor(center)"
                    else if (trait == "T2_Diabetes")
                        "Age + as.factor(Current_Smoker) + Systolic_Blood_Pressure + Total_Cholesterol + HDL_Cholesterol + as.factor(BP_medication) + as.factor(center)"
                    else
                        "Age + as.factor(sex) + as.factor(Current_Smoker) + as.factor(T2D) + Systolic_Blood_Pressure + Total_Cholesterol + HDL_Cholesterol + as.factor(BP_medication) + as.factor(center)"
                },
                PANEL = {
                    if(remove_sex)
                        "Age + as.factor(Current_Smoker) + 
                        as.factor(Alcohol_Intake) + Daily_Physical_Activity_minutes + Education +
                        as.factor(Daily_vegetable) + as.factor(FH_Type_2_Diabetes) + as.factor(T2D) + 
                        Body_Mass_Index + WHR + Waist_Circumference + Weight + Standing_Height + 
                        Systolic_Blood_Pressure + Total_Cholesterol + HDL_Cholesterol + LDL_Cholesterol +
                        Triglycerides + Glucose + HbA1c + Creatinine + Cystatine_C + Urea + Urate + 
                        AST + ALT + AP + Albumin + C_reaktive_Protein + Erythrocyte_Count + Leukocyte_Count +
                        Platelet_Count + Hemoglobin + Hematocrit + MCH + MCV + MCHC + as.factor(BP_medication) + as.factor(center)"
                    else if (trait == "T2_Diabetes")
                        "Age + as.factor(Current_Smoker) + 
                        as.factor(Alcohol_Intake) + Daily_Physical_Activity_minutes + Education +
                        as.factor(Daily_vegetable) + as.factor(FH_Type_2_Diabetes) + 
                        Body_Mass_Index + WHR + Waist_Circumference + Weight + Standing_Height + 
                        Systolic_Blood_Pressure + Total_Cholesterol + HDL_Cholesterol + LDL_Cholesterol +
                        Triglycerides + Glucose + HbA1c + Creatinine + Cystatine_C + Urea + Urate + 
                        AST + ALT + AP + Albumin + C_reaktive_Protein + Erythrocyte_Count + Leukocyte_Count +
                        Platelet_Count + Hemoglobin + Hematocrit + MCH + MCV + MCHC + as.factor(BP_medication) + as.factor(center)"
                    else
                        "Age + as.factor(sex) + as.factor(Current_Smoker) + 
                        as.factor(Alcohol_Intake) + Daily_Physical_Activity_minutes + Education +
                        as.factor(Daily_vegetable) + as.factor(FH_Type_2_Diabetes) + as.factor(T2D) + 
                        Body_Mass_Index + WHR + Waist_Circumference + Weight + Standing_Height + 
                        Systolic_Blood_Pressure + Total_Cholesterol + HDL_Cholesterol + LDL_Cholesterol +
                        Triglycerides + Glucose + HbA1c + Creatinine + Cystatine_C + Urea + Urate + 
                        AST + ALT + AP + Albumin + C_reaktive_Protein + Erythrocyte_Count + Leukocyte_Count +
                        Platelet_Count + Hemoglobin + Hematocrit + MCH + MCV + MCHC + as.factor(BP_medication) + as.factor(center)"
                },
                stop("Unknown covar.set.type")
            )
    Cs.cov <- c()
    Cs.m <- c()
    Cs.p <- c()
    Cs.mp <- c()
    for (rs in randomseeds){
        # read idsplit
        idsplit = fread(paste0("idsplit_overlap_",rs,".tab"))
        idsplit = idsplit %>% select(f.eid, !!sym(trait)) %>% rename(batch := !!sym(trait))
        rs.covar.df = covar.df %>% right_join(idsplit,by="f.eid") # get only overlap individuals

        mp.res = fread(paste0("mp_",trait,"_",covar.set.type,"_", rs, ".txt"))
        m.res = fread(paste0("/mp_meta_",trait, "_",covar.set.type,"_", rs, ".txt"))
        p.res = fread(paste0("/mp_prote_",trait,"_",covar.set.type,"_", rs, ".txt"))
        mp.res.df1 = rs.covar.df %>% right_join(mp.res %>% select(f.eid,pred),by="f.eid")
        m.res.df1 = rs.covar.df %>% right_join(m.res %>% select(f.eid,pred),by="f.eid")
        p.res.df1 = rs.covar.df %>% right_join(p.res %>% select(f.eid,pred),by="f.eid")
      
        # covariates only model
        cs = c()
        for (b in c(1:5)){
            train.covar.df = mp.res.df1 %>% filter(batch != b)
            test.covar.df = mp.res.df1 %>% filter(batch == b)
            fit1.formula = as.formula(paste0("Surv(event.age, train.covar.df[[trait]]) ~ ",covar.set))
            fit1 = coxph(fit1.formula, data = train.covar.df)
            risk = predict(fit1, newdata = test.covar.df, type = "risk")
            c = concordance.index(x = risk, surv.time = test.covar.df$event.age, surv.event = test.covar.df[[trait]],method = "noether")
            cs = c(cs,c$c.index)
        }
        Cs.cov = c(Cs.cov,mean(cs))
      
        # add metabolites only residuals
        fit2.formula = as.formula(paste0("Surv(event.age, m.res.df1[[trait]]) ~ ",covar.set, " + pred"))
        fit2 = coxph(fit2.formula, data = m.res.df1)
        Cs.m = c(Cs.m, summary(fit2)$concordance[1])
      
        # add proteomics only residuals
        fit3.formula = as.formula(paste0("Surv(event.age, p.res.df1[[trait]]) ~ ",covar.set, " + pred"))
        fit3 = coxph(fit3.formula, data = p.res.df1)
        Cs.p = c(Cs.p, summary(fit3)$concordance[1])
      
        # add metabolites + proteomics residuals
        fit4.formula = as.formula(paste0("Surv(event.age, mp.res.df1[[trait]]) ~ ",covar.set, " + pred"))
        fit4 = coxph(fit4.formula, data = mp.res.df1)
        Cs.mp = c(Cs.mp, summary(fit4)$concordance[1])
    }

    return(c(trait, c(Cs.cov), c(Cs.m), c(Cs.p), c(Cs.mp)))
}

result1 = NULL
result2 = NULL
result3 = NULL
kept_traits = c('Dementia', 'MACE',
        'T2_Diabetes', 'Liver_Disease', 'Renal_Disease', 'Atrial_Fibrillation',
        'Heart_Failure', 'CHD', 
        'PAD', 'Asthma', 'COPD', 'Skin_Cancer', 
        'Prostate_Cancer', 'Breast_Cancer', 
        'Fractures', 'Cataracts', 'Glaucoma')
for (t in kept_traits){
    result1 = rbind(result1,adj(comb.df,t,"age_sex"))
    result2 = rbind(result2,adj(comb.df,t, "ASCVD"))
    result3 = rbind(result3,adj(comb.df,t, "PANEL"))
    print(paste0(t," done"))
}

pd = function(result){
    result = as.data.frame(result)
    colnames(result) = c("Trait", "C.covar.fold1","C.covar.fold2",
                         "C.covar.fold3","C.covar.fold4","C.covar.fold5",
                         "C.meta.fold1","C.meta.fold2","C.meta.fold3",
                         "C.meta.fold4","C.meta.fold5",
                         "C.prote.fold1","C.prote.fold2","C.prote.fold3",
                         "C.prote.fold4","C.prote.fold5",
                         "C.mp.fold1","C.mp.fold2","C.mp.fold3",
                         "C.mp.fold4","C.mp.fold5")
    return(result)
}

result1 = pd(result1)
result2 = pd(result2)
result3 = pd(result3)
longformat = function(df,covar){
    df.long = df %>%
        pivot_longer(
            cols = starts_with("C."),
            names_to = c("model", "fold"),
            names_pattern = "C\\.(.*)\\.fold(.*)",
            values_to = "value"
        ) %>% mutate(model = case_when(
            model == "covar" ~ "Baseline covariates",
            model == "meta" ~ "Metabolomics-only",
            model == "prote" ~ "Proteomics-only",
            model == "mp" ~ "Combined"
        ))
    write.table(df.long, paste0("all_results_",covar,"_5rs_tables.txt"), quote=F, row.names=F, sep="\t")
}
longformat(result1, "age_sex")
longformat(result2, "ASCVD")
longformat(result3, "PANEL")

df1 = result1 %>%
    mutate(across(contains(".fold"), as.numeric)) %>%
    mutate(
        `mean.Age+Sex` = rowMeans(select(., starts_with("C.covar.fold")), na.rm = TRUE),
        `sd.Age+Sex`   = apply(select(., starts_with("C.covar.fold")), 1, sd, na.rm = TRUE),
        `min.Age+Sex`  = apply(select(., starts_with("C.covar.fold")), 1, min, na.rm = TRUE),
        `max.Age+Sex`  = apply(select(., starts_with("C.covar.fold")), 1, max, na.rm = TRUE),
        
        `mean.Meta+Age+Sex` = rowMeans(select(., starts_with("C.meta.fold")), na.rm = TRUE),
        `sd.Meta+Age+Sex`   = apply(select(., starts_with("C.meta.fold")), 1, sd, na.rm = TRUE),
        `min.Meta+Age+Sex`  = apply(select(., starts_with("C.meta.fold")), 1, min, na.rm = TRUE),
        `max.Meta+Age+Sex`  = apply(select(., starts_with("C.meta.fold")), 1, max, na.rm = TRUE),
        
        `mean.Proteo+Age+Sex` = rowMeans(select(., starts_with("C.prote.fold")), na.rm = TRUE),
        `sd.Proteo+Age+Sex`   = apply(select(., starts_with("C.prote.fold")), 1, sd, na.rm = TRUE),
        `min.Proteo+Age+Sex`  = apply(select(., starts_with("C.prote.fold")), 1, min, na.rm = TRUE),
        `max.Proteo+Age+Sex`  = apply(select(., starts_with("C.prote.fold")), 1, max, na.rm = TRUE),
        
        `mean.Meta+Proteo+Age+Sex` = rowMeans(select(., starts_with("C.mp.fold")), na.rm = TRUE),
        `sd.Meta+Proteo+Age+Sex`   = apply(select(., starts_with("C.mp.fold")), 1, sd, na.rm = TRUE),
        `min.Meta+Proteo+Age+Sex`  = apply(select(., starts_with("C.mp.fold")), 1, min, na.rm = TRUE),
        `max.Meta+Proteo+Age+Sex`  = apply(select(., starts_with("C.mp.fold")), 1, max, na.rm = TRUE)
    ) %>%
    select(Trait, starts_with("mean."), starts_with("sd."), starts_with("min."), starts_with("max."))

df2 = result2 %>%
    mutate(across(contains(".fold"), as.numeric)) %>%
    mutate(
        `mean.ASCVD` = rowMeans(select(., starts_with("C.covar.fold")), na.rm = TRUE),
        `sd.ASCVD`   = apply(select(., starts_with("C.covar.fold")), 1, sd, na.rm = TRUE),
        `min.ASCVD`  = apply(select(., starts_with("C.covar.fold")), 1, min, na.rm = TRUE),
        `max.ASCVD`  = apply(select(., starts_with("C.covar.fold")), 1, max, na.rm = TRUE),
        
        `mean.Meta+ASCVD` = rowMeans(select(., starts_with("C.meta.fold")), na.rm = TRUE),
        `sd.Meta+ASCVD`   = apply(select(., starts_with("C.meta.fold")), 1, sd, na.rm = TRUE),
        `min.Meta+ASCVD`  = apply(select(., starts_with("C.meta.fold")), 1, min, na.rm = TRUE),
        `max.Meta+ASCVD`  = apply(select(., starts_with("C.meta.fold")), 1, max, na.rm = TRUE),
        
        `mean.Proteo+ASCVD` = rowMeans(select(., starts_with("C.prote.fold")), na.rm = TRUE),
        `sd.Proteo+ASCVD`   = apply(select(., starts_with("C.prote.fold")), 1, sd, na.rm = TRUE),
        `min.Proteo+ASCVD`  = apply(select(., starts_with("C.prote.fold")), 1, min, na.rm = TRUE),
        `max.Proteo+ASCVD`  = apply(select(., starts_with("C.prote.fold")), 1, max, na.rm = TRUE),
        
        `mean.Meta+Proteo+ASCVD` = rowMeans(select(., starts_with("C.mp.fold")), na.rm = TRUE),
        `sd.Meta+Proteo+ASCVD`   = apply(select(., starts_with("C.mp.fold")), 1, sd, na.rm = TRUE),
        `min.Meta+Proteo+ASCVD`  = apply(select(., starts_with("C.mp.fold")), 1, min, na.rm = TRUE),
        `max.Meta+Proteo+ASCVD`  = apply(select(., starts_with("C.mp.fold")), 1, max, na.rm = TRUE)
    ) %>%
    select(Trait, starts_with("mean."), starts_with("sd."), starts_with("min."), starts_with("max."))

df3 = result3 %>%
    mutate(across(contains(".fold"), as.numeric)) %>%
    mutate(
        `mean.PANEL` = rowMeans(select(., starts_with("C.covar.fold")), na.rm = TRUE),
        `sd.PANEL`   = apply(select(., starts_with("C.covar.fold")), 1, sd, na.rm = TRUE),
        `min.PANEL`  = apply(select(., starts_with("C.covar.fold")), 1, min, na.rm = TRUE),
        `max.PANEL`  = apply(select(., starts_with("C.covar.fold")), 1, max, na.rm = TRUE),
        
        `mean.Meta+PANEL` = rowMeans(select(., starts_with("C.meta.fold")), na.rm = TRUE),
        `sd.Meta+PANEL`   = apply(select(., starts_with("C.meta.fold")), 1, sd, na.rm = TRUE),
        `min.Meta+PANEL`  = apply(select(., starts_with("C.meta.fold")), 1, min, na.rm = TRUE),
        `max.Meta+PANEL`  = apply(select(., starts_with("C.meta.fold")), 1, max, na.rm = TRUE),
        
        `mean.Proteo+PANEL` = rowMeans(select(., starts_with("C.prote.fold")), na.rm = TRUE),
        `sd.Proteo+PANEL`   = apply(select(., starts_with("C.prote.fold")), 1, sd, na.rm = TRUE),
        `min.Proteo+PANEL`  = apply(select(., starts_with("C.prote.fold")), 1, min, na.rm = TRUE),
        `max.Proteo+PANEL`  = apply(select(., starts_with("C.prote.fold")), 1, max, na.rm = TRUE),
        
        `mean.Meta+Proteo+PANEL` = rowMeans(select(., starts_with("C.mp.fold")), na.rm = TRUE),
        `sd.Meta+Proteo+PANEL`   = apply(select(., starts_with("C.mp.fold")), 1, sd, na.rm = TRUE),
        `min.Meta+Proteo+PANEL`  = apply(select(., starts_with("C.mp.fold")), 1, min, na.rm = TRUE),
        `max.Meta+Proteo+PANEL`  = apply(select(., starts_with("C.mp.fold")), 1, max, na.rm = TRUE)
    ) %>%
    select(Trait, starts_with("mean."), starts_with("sd."), starts_with("min."), starts_with("max."))

combined = df1 %>% left_join(df2, by = "Trait") %>% left_join(df3, by = "Trait")

plot.list = list()
kept_traits_ordered = c('MACE','CHD', 'Heart_Failure', 'Atrial_Fibrillation','PAD','T2_Diabetes', 
        'Prostate_Cancer', 'Breast_Cancer', 'Skin_Cancer', 
        'Dementia', 'Cataracts', 'Glaucoma','Asthma', 'COPD', 'Liver_Disease', 'Renal_Disease','Fractures')
plot.list = list()
i=1
for (t in kept_traits_ordered){
    tmp = combined %>%
        filter(Trait == t) %>%
        pivot_longer(cols = -Trait,
                    names_to = c("metric", "covariates sets"),
                    names_pattern = "(.*)\\.(.*)",
                    values_to = "value") %>%
        pivot_wider(
                names_from = metric,
                values_from = value
        ) %>%
        mutate(
            mean = as.numeric(mean),
            min = as.numeric(min),
            max = as.numeric(max),
            sd = as.numeric(sd)
            ) %>%
        mutate(`covariates sets` = factor(`covariates sets`,
                levels = c("Age+Sex", "Meta+Age+Sex", "Proteo+Age+Sex", "Meta+Proteo+Age+Sex",
                "ASCVD", "Meta+ASCVD", "Proteo+ASCVD", "Meta+Proteo+ASCVD",
                "PANEL", "Meta+PANEL", "Proteo+PANEL", "Meta+Proteo+PANEL")),

        color = c(
        "blue", "green", "red", "purple",
        "forestgreen", "pink", "brown", "cyan",
        "magenta", "darkorange", "darkgray", "black"
        ))

    p = ggplot(data = tmp,aes(x = `covariates sets`,y = mean)) +
            geom_point(aes(color = color), size = 0.5) +
            geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = color), width = 0.4,linewidth=0.5) +
            scale_color_identity() +
            scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
            labs(title = t, x=NULL,y=NULL) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
                    axis.text.y = element_text(size = 14),
                    plot.title = element_text(hjust = 0.5))
        plot.list[[i]] = p
        i = i+1
}

combine.pdf = ggarrange(plotlist = plot.list, ncol = 5, nrow = 4)

ggsave(paste0("all_results_5rs_phase2.pdf"),combine.pdf,width=16,height=20)
