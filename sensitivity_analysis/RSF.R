# fit random survival forest
library(ranger)
library(dplyr)
library(data.table)
library(survival)
library(glmnet)
library(readr)
library(stringr)

meta_keep <- fread("meta_keep.csv")
colnames(meta_keep) <- gsub("[ -]", "_", colnames(meta_keep))
print("\nFinish reading metabolites")

prote = fread("proteomics_full.tab")
print("\nFinish reading proteomics")

covar.imp = fread("covariates_imp_new.tab")
tte = fread("ICD10_time_to_event.tab")
comb.df = tte %>% right_join(covar.imp,by="f.eid")
comb.df$baseline_visit = as.Date(comb.df$baseline_visit)
comb.df$DOB = as.Date(comb.df$DOB)
print("\nFinish reading covariates and time-to-event")


fit.RSF = function(covar.df, trait, covar.set.type, rs){
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
                        Triglycerides + Glucose_covar + HbA1c + Creatinine_covar + Cystatine_C + Urea + Urate + 
                        AST + ALT + AP + Albumin_covar + C_reaktive_Protein + Erythrocyte_Count + Leukocyte_Count +
                        Platelet_Count + Hemoglobin + Hematocrit + MCH + MCV + MCHC + as.factor(BP_medication) + as.factor(center)"
                    else
                        "Age + as.factor(sex) + as.factor(Current_Smoker) + 
                        as.factor(Alcohol_Intake) + Daily_Physical_Activity_minutes + Education +
                        as.factor(Daily_vegetable) + as.factor(FH_Type_2_Diabetes) + as.factor(T2D) + 
                        Body_Mass_Index + WHR + Waist_Circumference + Weight + Standing_Height + 
                        Systolic_Blood_Pressure + Total_Cholesterol + HDL_Cholesterol + LDL_Cholesterol +
                        Triglycerides + Glucose_covar + HbA1c + Creatinine_covar + Cystatine_C + Urea + Urate + 
                        AST + ALT + AP + Albumin_covar + C_reaktive_Protein + Erythrocyte_Count + Leukocyte_Count +
                        Platelet_Count + Hemoglobin + Hematocrit + MCH + MCV + MCHC + as.factor(BP_medication) + as.factor(center)"
                },
                stop("Unknown covar.set.type")
            )
    if (covar.set.type == "PANEL"){
        selected_covar_df = covar.df %>% rename(Glucose_covar = Glucose, Albumin_covar = Albumin, Creatinine_covar = Creatinine) %>% select(f.eid, event.age, !!sym(trait), all_of(all.vars(as.formula(paste("~", covar.set)))))
    } else {
        selected_covar_df = covar.df %>% select(f.eid, event.age, !!sym(trait), all_of(all.vars(as.formula(paste("~", covar.set)))))
    }
    results_list <- list()

    print(paste0("Processing ", rs))
    idsplit = fread(paste0("idsplit_overlap_",rs,".tab")) %>% select(f.eid, !!sym(trait)) %>% rename(batch := !!sym(trait))
    all_df = idsplit %>% left_join(selected_covar_df, by="f.eid") %>% left_join(meta_keep, by="f.eid") %>% left_join(prote, by="f.eid") %>% na.omit()
    print(dim(all_df))
    y_surv = Surv(all_df$event.age, all_df[[trait]])
    X_covar <- model.matrix(as.formula(paste0("~", covar.set)), data = all_df)[,-1]
    colnames(X_covar) <- make.names(colnames(X_covar))
    meta_cols <- setdiff(colnames(meta_keep), "f.eid")
    prote_cols <- setdiff(colnames(prote), "f.eid")
    X_meta <- as.matrix(all_df[, ..meta_cols])
    colnames(X_meta) = make.names(colnames(X_meta))
    X_prote <- as.matrix(all_df[, ..prote_cols])
    colnames(X_prote) = make.names(colnames(X_prote))
    X_covar_meta = cbind(X_covar, X_meta)
    X_covar_prote = cbind(X_covar, X_prote)
    X_all <- cbind(X_covar, X_meta, X_prote)
    pred_covar = numeric(nrow(all_df))
    pred_meta = numeric(nrow(all_df))
    pred_prote = numeric(nrow(all_df))
    pred_combined = numeric(nrow(all_df))

    for (i in c(1:5)){
        print(paste0("  Fold ", i))
        train_idx = which(all_df$batch != i)
        test_idx = which(all_df$batch == i)
        y_train <- y_surv[train_idx]

        # baseline
        covar_train_df <- as.data.frame(X_covar[train_idx, ])
        covar_train_df$time <- y_train[, 1]
        covar_train_df$status <- y_train[, 2]
        fit_covar <- ranger(formula = Surv(time, status) ~ ., data = covar_train_df, num.trees = 1000, mtry = floor(sqrt(ncol(X_covar))), splitrule = "logrank")
        pred_covar[test_idx] <- rowSums(predict(fit_covar, data = as.data.frame(X_covar[test_idx, ]))$chf)

        # covar+meta 
        covar_meta_train_df <- as.data.frame(X_covar_meta[train_idx, ])
        covar_meta_train_df$time <- y_train[, 1]
        covar_meta_train_df$status <- y_train[, 2]
        covar.meta.fit = ranger(formula = Surv(time, status) ~ ., data = covar_meta_train_df, num.trees = 500, mtry = floor(sqrt(ncol(X_covar_meta))), splitrule = "logrank")
        pred_meta[test_idx] <- rowSums(predict(covar.meta.fit, data = as.data.frame(X_covar_meta[test_idx, ]))$chf)
        
        # covar+prote 
        covar_prote_train_df <- as.data.frame(X_covar_prote[train_idx, ])
        covar_prote_train_df$time <- y_train[, 1]
        covar_prote_train_df$status <- y_train[, 2]
        covar.prote.fit = ranger(formula = Surv(time, status) ~ ., data = covar_prote_train_df, num.trees = 500, mtry = floor(sqrt(ncol(X_covar_prote))), splitrule = "logrank")
        pred_prote[test_idx] <- rowSums(predict(covar.prote.fit, data = as.data.frame(X_covar_prote[test_idx, ]))$chf)

        # combined 
        all_train_df <- as.data.frame(X_all[train_idx, ])
        all_train_df$time <- y_train[, 1]
        all_train_df$status <- y_train[, 2]
        all.fit = ranger(formula = Surv(time, status) ~ ., data = all_train_df, num.trees = 500, mtry = floor(sqrt(ncol(X_all))), splitrule = "logrank")
        pred_combined[test_idx] <- rowSums(predict(all.fit, data = as.data.frame(X_all[test_idx, ]))$chf)

    }
    res_df <- data.frame(
        seed = rs,
        f.eid = all_df$f.eid,
        preds_covar = pred_covar,
        preds_meta = pred_meta,
        preds_prote = pred_prote,
        preds_combined = pred_combined
    )
    return(res_df)
}

t = commandArgs(trailingOnly=TRUE)[1]
covar_set_type = commandArgs(trailingOnly=TRUE)[2]
rs = commandArgs(trailingOnly=TRUE)[3]
result = NULL

print(paste0("Fitting trait: ", t))
result <- fit.RSF(comb.df, t, covar_set_type, rs)

fwrite(result, paste0("RSF_results/",t,"_",covar_set_type,"_",rs, "_pred.txt"),sep="\t")
