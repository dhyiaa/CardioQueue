#!/usr/bin/env Rscript

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

root <- "datasets/cardioboost/public_dataset/official_cardioBoost_model_snapshot"
rawroot <- "datasets/cardioboost/public_dataset/raw_public_cardioBoost_repo"
outdir <- "results/model_performance/cardioboost_matched_benchmark/tables"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write_prediction <- function(panel) {
  if (panel == "cm") {
    rm(list = setdiff(ls(envir = .GlobalEnv), c("root", "rawroot", "outdir", "write_prediction", "panel")))
    load(file.path(rawroot, "data/cardiomyopathy/cm_all_rare_mutation.RData"), envir = .GlobalEnv)
    cm_all_rare$pathogenic <- as.factor(cm_all_rare$pathogenic)
    load(file.path(root, "data/cardiomyopathy/preprocess.RData"), envir = .GlobalEnv)
    source(file.path(root, "script/src/preprocess_test.R"))
    source(file.path(root, "script/src/predict.R"))
    load(file.path(root, "data/cardiomyopathy/ml/train_ada.RData"), envir = .GlobalEnv)
    pred <- predict.cardioboost(train_ada, cm_all_rare)
    panel_name <- "CM"
    output <- file.path(outdir, "official_cardioboost_cm_all_rare_predictions.tsv")
  } else if (panel == "arm") {
    rm(list = setdiff(ls(envir = .GlobalEnv), c("root", "rawroot", "outdir", "write_prediction", "panel")))
    load(file.path(rawroot, "data/arrhythmia/arm_all_rare_mutation.RData"), envir = .GlobalEnv)
    arm_all_rare$pathogenic <- as.factor(arm_all_rare$pathogenic)
    load(file.path(root, "data/arrhythmia/preprocess.RData"), envir = .GlobalEnv)
    source(file.path(root, "script/src/preprocess_test.R"))
    source(file.path(root, "script/src/predict.R"))
    load(file.path(root, "data/arrhythmia/ml/train_ada.RData"), envir = .GlobalEnv)
    pred <- predict.cardioboost(train_ada, arm_all_rare)
    panel_name <- "ARM"
    output <- file.path(outdir, "official_cardioboost_arm_all_rare_predictions.tsv")
  } else {
    stop("Unknown panel: ", panel)
  }

  pred$cardioboost_panel <- panel_name
  pred$variant_id <- paste(pred$CHROM, pred$POS, pred$REF, pred$ALT, sep = "-")
  keep <- intersect(
    c(
      "variant_id", "cardioboost_panel", "CHROM", "POS", "REF", "ALT", "gene",
      "HGVSc", "HGVSp", "AF_Adj", "gnomAD_AF", "MCAP", "REVEL",
      "pathogenicity", "pathogenic", "key"
    ),
    colnames(pred)
  )
  write.table(pred[, keep], output, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  cat(panel_name, "rows:", nrow(pred), "output:", output, "\n")
}

write_prediction("cm")
write_prediction("arm")

write_holdout_prediction <- function(panel) {
  if (panel == "cm") {
    rm(list = setdiff(ls(envir = .GlobalEnv), c("root", "rawroot", "outdir", "write_prediction", "write_holdout_prediction", "panel")))
    load(file.path(rawroot, "data/cardiomyopathy/cm_holdout_test.RData"), envir = .GlobalEnv)
    test$pathogenic <- as.factor(test$pathogenic)
    load(file.path(root, "data/cardiomyopathy/preprocess.RData"), envir = .GlobalEnv)
    source(file.path(root, "script/src/preprocess_test.R"))
    source(file.path(root, "script/src/predict.R"))
    load(file.path(root, "data/cardiomyopathy/ml/train_ada.RData"), envir = .GlobalEnv)
    pred <- predict.cardioboost(train_ada, test)
    panel_name <- "CM"
    output <- file.path(outdir, "official_cardioboost_cm_released_holdout_predictions.tsv")
  } else if (panel == "arm") {
    rm(list = setdiff(ls(envir = .GlobalEnv), c("root", "rawroot", "outdir", "write_prediction", "write_holdout_prediction", "panel")))
    load(file.path(rawroot, "data/arrhythmia/arm_holdout_test.RData"), envir = .GlobalEnv)
    test$pathogenic <- as.factor(test$pathogenic)
    load(file.path(root, "data/arrhythmia/preprocess.RData"), envir = .GlobalEnv)
    source(file.path(root, "script/src/preprocess_test.R"))
    source(file.path(root, "script/src/predict.R"))
    load(file.path(root, "data/arrhythmia/ml/train_ada.RData"), envir = .GlobalEnv)
    pred <- predict.cardioboost(train_ada, test)
    panel_name <- "ARM"
    output <- file.path(outdir, "official_cardioboost_arm_released_holdout_predictions.tsv")
  } else {
    stop("Unknown panel: ", panel)
  }

  pred$cardioboost_panel <- panel_name
  pred$variant_id <- paste(pred$CHROM, pred$POS, pred$REF, pred$ALT, sep = "-")
  keep <- intersect(
    c(
      "variant_id", "cardioboost_panel", "CHROM", "POS", "REF", "ALT", "gene",
      "HGVSc", "HGVSp", "pathogenic", "pathogenicity", "key"
    ),
    colnames(pred)
  )
  write.table(pred[, keep], output, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  cat(panel_name, "released holdout rows:", nrow(pred), "output:", output, "\n")
}

write_holdout_prediction("cm")
write_holdout_prediction("arm")
