# 12_nested_cv_biomarker.R
# Nested cross-validation for the 9-gene exploratory biomarker panel
# (AD vs Old), using regularized logistic regression (glmnet).
#
# Outer loop (5-fold): honest performance estimate — each fold holds out
#   test samples never used in that round's model fitting or tuning.
# Inner loop (via cv.glmnet, 5-fold): selects the regularization strength
#   (lambda) within each outer training set only.
#
# IMPORTANT CAVEAT (see README/METHODS Limitations section): this nested
# CV only protects against leakage in lambda tuning. It does NOT protect
# against the fact that the 9 candidate genes were themselves selected
# using this same dataset (see 10_biomarker_candidates.R) — a separate,
# unaddressed source of circularity. Resulting AUCs should not be read as
# validated biomarker performance.

library(dplyr)
library(glmnet)
library(pROC)
library(org.Hs.eg.db)

## --- 1. Load data (same objects as 10_biomarker_candidates.R) -----------
normalized_data <- readRDS("data/normalized_data.rds")
col.data <- readRDS("data/col_data.rds")

# Same 9 candidate genes identified in 10_biomarker_candidates.R
immune_genes <- c("IL4I1", "BCAT2", "KMO", "SDS", "GLUD1", "KYAT1", "SELE", "NLRP3", "ELANE")

# We need the gene_id (not just the symbol) to subset normalized_data,
# since its rownames are Ensembl IDs, not symbols
res_table_sig_ad_old <- readRDS("data/res_table_sig_ad_old.rds")
exclusive_ad_old_df <- readRDS("data/genes_exclusive_ad_vs_old.rds")

candidates <- res_table_sig_ad_old %>%
  filter(sub("\\..*", "", gene_id) %in% sub("\\..*", "", exclusive_ad_old_df$gene_id))
candidates$gene_id_clean <- sub("\\..*", "", candidates$gene_id)
candidates$gene_symbol <- mapIds(org.Hs.eg.db, keys = candidates$gene_id_clean,
                                 keytype = "ENSEMBL", column = "SYMBOL")
candidates <- candidates %>% filter(gene_symbol %in% immune_genes)

## --- 2. Build predictor matrix: samples x genes, Old/AD only -----------
gene_counts <- normalized_data[candidates$gene_id, ]
rownames(gene_counts) <- candidates$gene_symbol

model_data <- as.data.frame(t(gene_counts))
model_data$condition <- col.data$condition[match(rownames(model_data), rownames(col.data))]
model_data <- model_data %>% filter(condition %in% c("Old", "AD"))
model_data$condition <- factor(model_data$condition, levels = c("Old", "AD"))

table(model_data$condition)  # should be 10 Old, 12 AD

## --- 3. Nested CV function ----------------------------------------------

library(caret)

set.seed(42)

# 5 folds estratificados: caret respeta la proporción Old/AD en cada uno
outer_folds <- createFolds(model_data$condition, k = 5, list = TRUE, returnTrain = FALSE)

x_all <- as.matrix(model_data[, immune_genes])
y_all <- model_data$condition

fold_aucs_ridge <- numeric(5)

for (i in seq_along(outer_folds)) {
  test_idx  <- outer_folds[[i]]
  train_idx <- setdiff(seq_len(nrow(model_data)), test_idx)
  
  x_train <- x_all[train_idx, ]
  y_train <- y_all[train_idx]
  x_test  <- x_all[test_idx, ]
  y_test  <- y_all[test_idx]
  
  # Inner CV: cv.glmnet elige el lambda que mejor generaliza DENTRO del train
  cv_fit <- cv.glmnet(x_train, y_train, family = "binomial",
                      alpha = 0, nfolds = 5, type.measure = "deviance")
  
  # Predecimos el fold de test externo, que el modelo nunca vio
  pred_prob <- predict(cv_fit, newx = x_test, s = "lambda.min", type = "response")
  
  roc_obj <- roc(y_test, as.numeric(pred_prob), quiet = TRUE)
  fold_aucs_ridge[i] <- as.numeric(auc(roc_obj))
}

fold_aucs_ridge
mean(fold_aucs_ridge)


run_nested_cv <- function(alpha_value, seed = 42) {
  set.seed(seed)
  outer_folds <- createFolds(model_data$condition, k = 5, list = TRUE, returnTrain = FALSE)
  
  fold_aucs <- numeric(5)
  coef_list <- vector("list", 5)  # coeficientes del modelo final de cada fold
  
  for (i in seq_along(outer_folds)) {
    test_idx  <- outer_folds[[i]]
    train_idx <- setdiff(seq_len(nrow(model_data)), test_idx)
    
    x_train <- x_all[train_idx, ]
    y_train <- y_all[train_idx]
    x_test  <- x_all[test_idx, ]
    y_test  <- y_all[test_idx]
    
    cv_fit <- cv.glmnet(x_train, y_train, family = "binomial",
                        alpha = alpha_value, nfolds = 5, type.measure = "deviance")
    
    pred_prob <- predict(cv_fit, newx = x_test, s = "lambda.min", type = "response")
    roc_obj <- roc(y_test, as.numeric(pred_prob), quiet = TRUE)
    fold_aucs[i] <- as.numeric(auc(roc_obj))
    
    coef_list[[i]] <- as.matrix(coef(cv_fit, s = "lambda.min"))
  }
  
  list(aucs = fold_aucs, coefs = coef_list)
}


## --- 4. Run for Ridge and Lasso ------------------------------------------
ridge_result <- run_nested_cv(alpha_value = 0)
lasso_result <- run_nested_cv(alpha_value = 1)

# Note: "dangerous ground" warnings (class with <8 obs in some inner fold)
# are expected given the small sample size — documented as a limitation,
# not a bug.

## --- 5. Compare results ---------------------------------------------------
auc_comparison <- data.frame(
  fold = 1:5,
  ridge_auc = ridge_result$aucs,
  lasso_auc = lasso_result$aucs
)
auc_comparison

cat("Ridge mean AUC:", mean(ridge_result$aucs), "\n")
cat("Lasso mean AUC:", mean(lasso_result$aucs), "\n")

# Gene retention across folds under Lasso (coefficient != 0 = retained)
lasso_coefs <- do.call(cbind, lapply(lasso_result$coefs, function(m) m[, 1]))
colnames(lasso_coefs) <- paste0("fold", 1:5)
round(lasso_coefs, 3)


# RESULT: both Ridge and Lasso reach AUC = 1.0 in all 5 outer folds.
# This reflects the circularity caveat above (gene selection used the
# full dataset), not validated biomarker performance — see Limitations.
# Lasso consistently zeroes IL4I1 (5/5 folds) and SELE (4/5 folds),
# suggesting redundancy within the panel, subject to the same caveat.

## --- 6. Save results ------------------------------------------------------
saveRDS(list(auc_comparison = auc_comparison, lasso_coefs = lasso_coefs),
        "data/nested_cv_results.rds")


