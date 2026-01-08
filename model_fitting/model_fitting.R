library(mixOmics)
library(caret)
library(dplyr)
library(data.table)
library(dplyr)
ncomp <- 10

get_SSY = function(object){
    if (!is.null(object$Y)) {
        Y <- object$Y
    } else if (!is.null(object$indY)) {
        Y <- object$X[[object$indY]]
    }

    if (is.factor(Y) || is.character(Y)) {
        Y_numeric <- as.numeric(as.factor(Y)) - 1 
    } else {
        Y_numeric <- as.matrix(Y)
    }
    n_comp = object$ncomp[1]
    block_names <- names(object$X)
    if (!is.null(object$indY)) {
        block_names <- block_names[-object$indY]
    }

    weight_list <- list()
    for (block_name in block_names) {
        T_block <- object$variates[[block_name]]
        X_block <- object$X[[block_name]]
        cor_Y <- cor(Y_numeric, T_block, use = "pairwise.complete.obs")^2
        if (ncol(as.matrix(Y_numeric)) > 1) {
            SS_Y <- apply(cor_Y, 2, sum)
        } else {
            SS_Y <- as.vector(cor_Y)
        }
        if (sum(SS_Y) == 0) SS_Y <- rep(1e-10, length(SS_Y))

        block_weight <- sum(SS_Y[1:n_comp]) / ncol(Y_numeric)
        weight_list[[block_name]] <- block_weight
    }
    return(weight_list)
}

# run combind model
runall <- function(data_path, id_path, idsplit_path, res_path, cov, trait, opt_dir, randomseed){
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

    preds <- data.frame(row.names = character(), pred = numeric())
    design = matrix(0.1, ncol = 2, nrow = 2, 
                  dimnames = list(c("metab", "prot"), c("metab", "prot")))
    diag(design) <- 0

    for (j in c(1:5)){
        test_idx <- rownames(folds[folds[[trait]] == j, ,drop = FALSE])
        train_idx <- rownames(folds[folds[[trait]] != j, ,drop = FALSE])
        meta_train <- meta[train_idx, ]
        meta_test <- meta[test_idx, ]
        prote_train <- prote[train_idx, ]
        prote_test <- prote[test_idx, ]
        train_res <- residue[train_idx, ,drop = FALSE]
        test_res <- residue[test_idx, ,drop = FALSE]
        X <- list(metab = as.matrix(meta_train), prot = as.matrix(prote_train))
        Y <- as.matrix(train_res)

        model <- block.spls(X, Y, ncomp = ncomp, mode = "regression", scale = F, design = design) # already scaled
        pred <- predict(model, newdata = list(metab = meta_test, prot = prote_test))
        weights = get_SSY(model)
        meta_weights = weights$metab / (weights$metab + weights$prot)
        prot_weights = weights$prot / (weights$metab + weights$prot)
        pred_metab_raw <- pred$predict[["metab"]][, , ncomp]
        pred_prot_raw  <- pred$predict[["prot"]][, , ncomp]
        pred_Y_weighted <- (meta_weights * pred_metab_raw) + (prot_weights * pred_prot_raw)
        preds[test_idx, "pred"] <- as.numeric(pred_Y_weighted)
    }
    preds$f.eid <- rownames(preds)
    preds <- preds %>% select(f.eid, everything())
    write.csv(preds, file.path(opt_dir, paste0("mp_", trait, "_", cov, "_", randomseed, ".txt")), row.names = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
runall(args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8])
