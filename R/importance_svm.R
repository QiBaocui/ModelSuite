# Saved radial-SVM importance, adapted from model.4.0.R get_varimp_svm().
# Fixed model, mean AUC drop after permutation of each original variable.

svm_auc <- function(y, probability) {
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  if (!n1 || !n0 || any(!is.finite(probability))) stop("Invalid AUC inputs.")
  (sum(rank(probability, ties.method = "average")[y == 1L]) -
    n1 * (n1 + 1) / 2) / (n1 * n0)
}

format_svm_varimp <- function(scores) {
  maximum <- max(abs(scores))
  positive <- pmax(scores, 0)
  positive_total <- sum(positive)
  result <- data.frame(variable = names(scores), relative_importance = unname(scores),
    scaled_importance = if (maximum > 0) unname(scores) / maximum else 0,
    percentage = if (positive_total > 0) unname(positive) / positive_total else 0,
    stringsAsFactors = FALSE)
  result <- result[order(-result$relative_importance, result$variable), , drop = FALSE]
  rownames(result) <- NULL
  result
}

extract_svm_varimp <- function(
  models_root, data_file,
  evaluation_dataset = c("temporal", "development"), n_perm = 100L,
  seed = 2026L, schemes = NULL, outdir = NULL, overwrite = FALSE
) {
  evaluation_dataset <- match.arg(evaluation_dataset)
  if (length(n_perm) != 1L || !is.finite(n_perm) || n_perm < 2L || n_perm != floor(n_perm))
    stop("n_perm must be an integer >= 2.")
  n_perm <- as.integer(n_perm)
  if (!requireNamespace("e1071", quietly = TRUE)) stop("Install e1071 first.")
  root <- models_root
  if (is.null(schemes)) {
    candidates <- list.dirs(root, recursive = FALSE, full.names = FALSE)
    schemes <- candidates[vapply(candidates, function(s)
      file.exists(file.path(root, s, "development/manifest.rds")), logical(1))]
  }
  if (!length(schemes) || anyDuplicated(schemes)) stop("No schemes or duplicate schemes.")
  manifests <- setNames(lapply(schemes, function(s)
    read_model_manifest(file.path(root, s, "development"))), schemes)
  if (is.null(outdir)) outdir <- file.path(root, "svm_varimp", evaluation_dataset)
  if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)) && !overwrite)
    stop("Output directory is not empty; choose another outdir or set overwrite=TRUE.")
  data <- read.delim(data_file, check.names = FALSE, stringsAsFactors = FALSE,
    colClasses = setNames("character", manifests[[1L]]$config$id_col))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  results <- list()
  for (scheme in schemes) {
    man <- manifests[[scheme]]; cfg <- man$config; prep <- man$preprocessor
    if (!"SVM" %in% cfg$models || !file.exists(man$model_files$SVM))
      stop("SVM model absent: ", scheme)
    model <- readRDS(man$model_files$SVM)
    # e1071 stores kernel as an integer (2 = radial), not the string 'radial'.
    if (!identical(as.integer(model$kernel), 2L)) stop("Expected radial SVM: ", scheme)
    stage <- if (evaluation_dataset == "temporal") "temporal" else "development/apparent"
    saved <- read.delim(file.path(root, scheme, stage, "SVM/predictions.tsv"),
      stringsAsFactors = FALSE, colClasses = c(ID = "character"))
    if (anyDuplicated(as.character(data[[cfg$id_col]]))) stop("Duplicate input IDs.")
    index <- match(saved$ID, as.character(data[[cfg$id_col]]))
    if (anyNA(index)) stop("Saved patients missing from input data: ", scheme)
    evaluation <- data[index, , drop = FALSE]
    y <- check_data(evaluation, cfg$id_col, cfg$outcome_col, cfg$predictors, cfg$positive, cfg$negative)
    if (!identical(as.integer(saved$y), y)) stop("Saved outcome mismatch: ", scheme)
    raw <- evaluation[, cfg$predictors, drop = FALSE]
    # Apply the original fitted imputation, dummy coding and scaling every time.
    predict_raw <- function(d) predict_engine("SVM", model, apply_preprocessor(prep, d))
    base_probability <- predict_raw(raw)
    if (max(abs(base_probability - saved$probability)) > 1e-7)
      stop("Baseline predictions differ from saved predictions: ", scheme)
    base_auc <- svm_auc(y, base_probability)
    set.seed(seed)
    scores <- setNames(vapply(cfg$predictors, function(variable) {
      mean(replicate(n_perm, {
        permuted <- raw
        # Jointly alter all encoded levels by permuting the original column.
        permuted[[variable]] <- raw[[variable]][sample.int(nrow(raw))]
        base_auc - svm_auc(y, predict_raw(permuted))
      }))
    }, numeric(1)), cfg$predictors)
    result <- format_svm_varimp(scores)
    folder <- file.path(outdir, scheme)
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    write.table(result, file.path(folder, paste0(scheme, "_SVM_varimp.tsv")),
      sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
    results[[scheme]] <- result
    message("Exported: ", scheme, "; N=", length(y), "; baseline AUC=",
      signif(base_auc, 6), "; variables=", nrow(result))
  }
  writeLines(c(
    "SVM permutation variable importance (saved final radial models; no refitting).",
    paste0("Evaluation dataset: ", evaluation_dataset, "; permutations per variable: ", n_perm, "; seed: ", seed, "."),
    "Each scheme folder contains <scheme>_SVM_varimp.tsv; four columns match the H2O table schema.",
    "variable: original predictor name; all categorical levels are permuted jointly.",
    "relative_importance: mean(baseline AUC - permuted AUC), using positive-class probabilities.",
    "AUC uses fixed positive/negative labels and direction, including when a permuted AUC is below 0.5.",
    "scaled_importance: relative_importance / max(abs(relative_importance)); zero when all scores are zero.",
    "percentage: max(relative_importance, 0) / sum(max(relative_importance, 0)), a 0-1 fraction; all zero when no score is positive.",
    "Negative raw/scaled scores are retained; negative scores receive percentage=0.",
    "Compared with model.4.0.R: 100 permutations by default (rather than 5), original-variable grouping, corrected scalar-ifelse percentage bug, and 0-1 fractions rather than 0-100.",
    "Temporal scores describe the locked model on the independent cohort; they are exploratory explanations and must not be used to retune the validated model.",
    "Development scores explain the final model on its own training patients, not OOF validation importance.",
    "Scores depend on the evaluation data and correlated predictors; they are not independent effects, causal effects or directly comparable with H2O native importance."
  ), file.path(outdir, "README.txt"))
  message("Complete: ", normalizePath(outdir, winslash = "/"))
  invisible(results)
}
