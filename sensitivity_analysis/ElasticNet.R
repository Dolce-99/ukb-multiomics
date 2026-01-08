# get prediction from elastic net
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

randomseeds <- c("rs42", "rs123", "rs456", "rs1024", "rs2048")

get_best_lambda_bic <- function(fit, n_samples) {
    bic_values <- fit$nulldev * (1 - fit$dev.ratio) + fit$df * log(n_samples)
    best_idx <- which.min(bic_values)
  
    return(list(
        lambda = fit$lambda[best_idx],
        idx = best_idx,
        bic_min = bic_values[best_idx]
    ))
}

fit.elasticNet <- function(trait, covar.set.type) {
    results_list <- list()
    res_name = paste0(trait,"_surv_res")
    res = fread(paste0("surv_res_",covar.set.type,"_overlap.tab")) %>% select(f.eid, !!sym(res_name))
    
    for (rs in randomseeds) {
        print(paste0("Processing ", rs))
        idsplit = fread(paste0("idsplit_overlap_",rs,".tab")) %>% select(f.eid, !!sym(trait)) %>% rename(batch := !!sym(trait))
        all_df = res %>% left_join(idsplit, by="f.eid") %>% left_join(meta_keep, by="f.eid") %>% left_join(prote, by="f.eid") %>% na.omit()
        print(dim(all_df))
        y = all_df[[res_name]]
        meta_cols <- setdiff(colnames(meta_keep), "f.eid")
        prote_cols <- setdiff(colnames(prote), "f.eid")
        X_meta <- as.matrix(all_df[, ..meta_cols])
        X_prote <- as.matrix(all_df[, ..prote_cols])
        X_all <- cbind(X_meta, X_prote)
        pred_meta = numeric(nrow(all_df))
        pred_prote = numeric(nrow(all_df))
        pred_combined = numeric(nrow(all_df))

        for (i in c(1:5)){
            print(paste0("  Fold ", i))
            train_idx = which(all_df$batch != i)
            test_idx = which(all_df$batch == i)
            y_train <- y[train_idx]
            n_train <- length(y_train)

            # meta bic
            meta.bic.fit = glmnet(x=X_meta[train_idx, ], y=y_train, family="gaussian", alpha=0.5)
            bic.meta = get_best_lambda_bic(meta.bic.fit, n_train)
            pred_meta[test_idx] = predict(meta.bic.fit, newx = X_meta[test_idx, ], s = bic.meta$lambda, type = "response")

            # prote bic
            prote.bic.fit = glmnet(x=X_prote[train_idx, ], y=y_train, family="gaussian", alpha=0.5)
            bic.prote = get_best_lambda_bic(prote.bic.fit, n_train)
            pred_prote[test_idx] = predict(prote.bic.fit, newx = X_prote[test_idx, ], s = bic.prote$lambda, type = "response")

            # combined bic
            combined.bic.fit = glmnet(x=X_all[train_idx, ], y=y_train, family="gaussian", alpha=0.5)
            bic.all = get_best_lambda_bic(combined.bic.fit, n_train)
            pred_combined[test_idx] = predict(combined.bic.fit, newx = X_all[test_idx, ], s = bic.all$lambda, type = "response")
        }
        res_df <- data.frame(
            seed = rs,
            f.eid = all_df$f.eid,
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
res <- fit.elasticNet(t, covar_set_type)
result <- rbind(result, res)

fwrite(result, paste0("elastic.net_results/",t,"_",covar_set_type,"_pred.txt"),sep="\t")




