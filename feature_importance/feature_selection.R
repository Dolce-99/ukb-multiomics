library(mixOmics)
library(caret)
library(dplyr)
ncomp <- 10

args <- commandArgs(trailingOnly = TRUE)

fit.model = function(data_path, omics, id_path, idsplit_path, res_path, cov, trait){
    data <- na.omit(read.table(paste0(data_path,"/",omics,".tab"), header = TRUE, sep = "\t", row.names = 1))
    ids <- read.table(id_path, header = TRUE, sep = "\t")$f.eid
  
    residue <- na.omit(read.table(res_path, header = TRUE, row.names = 1, sep = "\t")[, paste0(trait, "_surv_res"), drop = FALSE])
    folds <- read.table(idsplit_path, header = TRUE, row.names = 1, sep = "\t")[, trait, drop = FALSE]
    folds[[trait]] <- as.integer(folds[[trait]])
  
    common_id <- Reduce(intersect, list(ids, rownames(folds), rownames(residue)))
    data = data[common_id, ,drop = FALSE]
    folds <- folds[common_id, ,drop = FALSE]
    residue <- residue[common_id, ,drop = FALSE]

    j=1
    test_idx <- rownames(folds[folds[[trait]] == j, ,drop = FALSE])
    train_idx <- rownames(folds[folds[[trait]] != j, ,drop = FALSE])
    data_train <- data[train_idx, ]
    train_res <- residue[train_idx, ,drop = FALSE]
    X <- list(data = data_train)
    Y <- as.matrix(train_res)
    
    model <- block.spls(X, Y, ncomp = ncomp, mode = "regression", scale = F) # already scaled
    loadings <- model$loadings$data
    AVE <- unlist(model$AVE$AVE_X$data, use.names = FALSE)
    num_feat = ncol(data_train)
    weighted_W <- (sqrt(num_feat) * abs(loadings) %*% AVE) / sum(AVE)
    res_df = data.frame(feature = colnames(data), Omic = omics, Imp_raw = as.numeric(weighted_W))
    return(res_df)
}

meta_feat_df = fit.model(args[1], "meta", args[2], args[3], args[4], args[5], args[6])
prote_feat_df = fit.model(args[1], "prote", args[2], args[3], args[4], args[5], args[6])

runall <- function(data_path, id_path, idsplit_path, res_path, cov, trait, opt_dir, meta_feat_df, prote_feat_df){
    meta = na.omit(read.table(paste0(data_path,"/meta.tab"), header = TRUE, sep = "\t", row.names = 1))
    prote <- na.omit(read.table(paste0(data_path,"/prote.tab"), header = TRUE, sep = "\t", row.names = 1))
    ids <- read.table(id_path, header = TRUE, sep = "\t")$f.eid
  
    residue <- na.omit(read.table(res_path, header = TRUE, row.names = 1, sep = "\t")[, paste0(trait, "_surv_res"), drop = FALSE])
    folds <- read.table(idsplit_path, header = TRUE, row.names = 1, sep = "\t")[, trait, drop = FALSE]
    folds[[trait]] <- as.integer(folds[[trait]])
  
    common_id <- Reduce(intersect, list(ids, rownames(folds), rownames(residue)))
    meta = meta[common_id, ,drop = FALSE]
    prote = prote[common_id, ,drop = FALSE]
    folds <- folds[common_id, ,drop = FALSE]
    residue <- residue[common_id, ,drop = FALSE]

    j=1
    test_idx <- rownames(folds[folds[[trait]] == j, ,drop = FALSE])
    train_idx <- rownames(folds[folds[[trait]] != j, ,drop = FALSE])
    meta_train <- meta[train_idx, ]
    meta_test <- meta[test_idx, ]
    prote_train <- prote[train_idx, ]
    prote_test <- prote[test_idx, ]
    train_res <- residue[train_idx, ,drop = FALSE]
    test_res <- residue[test_idx, ,drop = FALSE]
    X <- list(meta = as.matrix(meta_train), prote = as.matrix(prote_train))
    design = matrix(0.1, ncol = 2, nrow = 2, 
                  dimnames = list(c("meta", "prote"), c("meta", "prote")))
    diag(design) <- 0
    Y <- as.matrix(train_res)
    
    model <- block.spls(X, Y, ncomp = ncomp, mode = "regression", scale = F, design = design) 
    
    final_results <- list()
    block_names <- names(X)
    imp_all = rbind(meta_feat_df, prote_feat_df)
    rownames(imp_all) = imp_all$feature

    for (block_name in block_names) {
        T_block <- model$variates[[block_name]]
        X_block <- X[[block_name]]
        cor_Y <- cor(Y, T_block, use = "pairwise.complete.obs")^2
        if (ncol(Y) > 1) { SS_Y_vec <- apply(cor_Y, 2, sum) } else { SS_Y_vec <- as.vector(cor_Y) }
        if (sum(SS_Y_vec) == 0) SS_Y_vec <- rep(1e-10, length(SS_Y_vec))
        block_weight <- sum(SS_Y_vec[1:ncomp]) / ncol(Y)
        common_rows <- intersect(rownames(X_block), rownames(Y))
        feat_cor <- cor(X_block[common_rows, ], Y[common_rows, ], use = "pairwise.complete.obs")
        features_in_block <- colnames(X_block)
        block_imp_data <- imp_all[features_in_block, ] %>% as.data.frame()
        imp_raw = block_imp_data$Imp_raw
        res_df <- data.frame(
            Feature = features_in_block,
            Block = block_name,
            Imp_Raw = imp_raw, # from single omics model
            Block_Weight = block_weight,     
            Correlation = as.numeric(feat_cor), 
            Imp_Weighted = imp_raw * block_weight,
            Final_Score = (imp_raw * block_weight) * sign(as.numeric(feat_cor)),  
            stringsAsFactors = FALSE
        )
        final_results[[block_name]] <- res_df
    }
    master_df <- do.call(rbind, final_results)
    master_df <- master_df[order(abs(master_df$Final_Score), decreasing = TRUE), ]
    write.table(master_df, file.path(opt_dir, paste0("mp_", trait, "_", cov, "_weighted_abs.txt")), row.names = TRUE,quote=FALSE, sep = "\t")
}

runall(args[1], args[2], args[3], args[4], args[5], args[6], args[7], meta_feat_df, prote_feat_df)

