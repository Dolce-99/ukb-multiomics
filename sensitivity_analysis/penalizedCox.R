library(dplyr)
library(data.table)
library(survival)
library(glmnet)
library(readr)
library(stringr)

meta <- fread("ukb_metabofilter_complete_covars.txt")
meta_field <- fread("metabolites.data.field.txt")
map_col_dict <- setNames(meta_field$Abbreviation, meta_field$`Field ID`)
update_col_names <- function(col_names) {
    sapply(col_names, function(x) {
        match_num <- str_extract(x, "\\d+")
    if (!is.na(match_num) && match_num %in% names(map_col_dict)) {
      return(map_col_dict[[match_num]])
    }
    return(x)
  })
}
colnames(meta) <- unname(update_col_names(colnames(meta)))
valid_metabolites <- map_col_dict[map_col_dict %in% colnames(meta)]
cols_to_keep <- c("n_eid", valid_metabolites)
meta_keep <- meta %>% select(any_of(cols_to_keep)) %>% rename(f.eid = n_eid)
print("\nFinish reading metabolites")

prote = fread("proteomics_full.tab")
print("\nFinish reading proteomics")

covar.imp = fread("covariates_imp_new.tab")
tte = fread("ICD10_time_to_event.tab")
comb.df = tte %>% right_join(covar.imp,by="f.eid")
comb.df$baseline_visit = as.Date(comb.df$baseline_visit)
comb.df$DOB = as.Date(comb.df$DOB)
print("\nFinish reading covariates and time-to-event")

randomseeds <- c("rs42", "rs123", "rs456", "rs1024", "rs2048")

get_best_lambda_bic <- function(fit, y_surv) {
    bic_values <- fit$nulldev * (1 - fit$dev.ratio) + fit$df * log(sum(y_surv[, 2]))
    best_idx <- which.min(bic_values)
  
    return(list(
        lambda = fit$lambda[best_idx],
        idx = best_idx,
        bic_min = bic_values[best_idx]
    ))
}

fit.penalized.cox <- function(covar.df, trait, covar.set.type) {
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
    selected_covar_df = covar.df %>% select(f.eid, event.age, !!sym(trait), all_of(all.vars(as.formula(paste("~", covar.set)))))
    results_list <- list()

    for (rs in randomseeds) {
        print(paste0("Processing ", rs))
        idsplit = fread(paste0("idsplit_overlap_",rs,".tab")) %>% select(f.eid, !!sym(trait)) %>% rename(batch := !!sym(trait))
        all_df = idsplit %>% left_join(selected_covar_df, by="f.eid") %>% left_join(meta_keep, by="f.eid") %>% left_join(prote, by="f.eid") %>% na.omit()
        print(dim(all_df))
        y_surv = Surv(all_df$event.age, all_df[[trait]])
        X_covar = model.matrix(as.formula(paste0("~", covar.set)), data=all_df)[,-1]
        meta_cols <- setdiff(colnames(meta_keep), "f.eid")
        prote_cols <- setdiff(colnames(prote), "f.eid")
        X_meta <- as.matrix(all_df[, ..meta_cols])
        X_prote <- as.matrix(all_df[, ..prote_cols])
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
            X_train <- X_all[train_idx, ]
            y_train <- y_surv[train_idx]
            X_test  <- X_all[test_idx, ]

            # baseline
            covar_train_df <- as.data.frame(X_covar[train_idx, ])
            covar_train_df$time <- y_train[, 1]
            covar_train_df$status <- y_train[, 2]
            fit_covar <- coxph(Surv(time, status) ~ ., data = covar_train_df)
            pred_covar[test_idx] <- predict(fit_covar, newdata = as.data.frame(X_covar[test_idx, ]))

            # covar+meta bic, lasso
            meta.bic.fit = glmnet(x=X_covar_meta[train_idx, ], y=y_train, family="cox", alpha=1, penalty.factor=c(rep(0, ncol(X_covar)), rep(1, ncol(X_meta))))
            bic.meta = get_best_lambda_bic(meta.bic.fit, y_train)
            pred_meta[test_idx] <- predict(meta.bic.fit, newx = X_covar_meta[test_idx, ], s = bic.meta$lambda, type = "link")

            # covar+prote bic, lasso
            prote.bic.fit = glmnet(x=X_covar_prote[train_idx, ], y=y_train, family="cox", alpha=1, penalty.factor=c(rep(0, ncol(X_covar)), rep(1, ncol(X_prote))))
            bic.prote = get_best_lambda_bic(prote.bic.fit, y_train)
            pred_prote[test_idx] <- predict(prote.bic.fit, newx = X_covar_prote[test_idx, ], s = bic.prote$lambda, type = "link")

            # combined bic, lasso
            combined.bic.fit = glmnet(x=X_all[train_idx, ], y=y_train, family="cox", alpha=1, penalty.factor=c(rep(0, ncol(X_covar)),rep(1, ncol(X_meta) + ncol(X_prote))))
            bic.all = get_best_lambda_bic(combined.bic.fit, y_train)
            pred_combined[test_idx] <- predict(combined.bic.fit, newx = X_test, s = bic.all$lambda, type = "link")
        }
        res_df <- data.frame(
            seed = rs,
            f.eid = all_df$f.eid,
            preds_covar = pred_covar,
            preds_meta = pred_meta,
            preds_prote = pred_prote,
            preds_combined = pred_combined
        )
        results_list[[rs]] <- res_df
    }
    return(bind_rows(results_list))
}


t = commandArgs(trailingOnly=TRUE)[1]
covar_set_type = commandArgs(trailingOnly=TRUE)[2]
result = NULL

print(paste0("Fitting trait: ", t))
res <- fit.penalized.cox(comb.df, t, covar_set_type)
result <- rbind(result, res)

fwrite(result, paste0("penalizedCox_results/",t,"_",covar_set_type,"_pred_LassoBIC.txt"),sep="\t")



