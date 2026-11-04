# Unified orchestration; stages reuse the same saved models and preprocessing.
run_pipeline <- function(cfg, stage = c("all", "check", "model", "importance", "plot")) {
  stage <- match.arg(stage)
  schemes <- names(cfg$schemes)
  models_root <- cfg$models_root
  importance_models <- cfg$importance_models
  if (stage %in% c("all", "check", "model")) {
    cohorts <- read_cohorts(cfg)
    message("Input checks passed: development=", nrow(cohorts$train), "; temporal=", nrow(cohorts$test),
      "; schemes=", length(schemes), "; models=", paste(cfg$models, collapse = ","))
    if (stage == "check") return(invisible(TRUE))
    # Validate all output destinations before starting the first expensive fit.
    for (scheme in schemes) for (part in c("development", "temporal")) {
      path <- file.path(models_root, scheme, part)
      if (dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE)))
        stop("Existing modeling output: ", path, "; use a new output_dir.")
    }
    for (scheme in schemes) {
      specification <- cfg$schemes[[scheme]]
      modelFun(train.R = cohorts$train, test.R = cohorts$test, outdir = models_root, feature = scheme,
        feature_type_list = cfg$models, outcome_col = cfg$outcome_col, id_col = cfg$id_col,
        positive = cfg$positive, predictors = specification$predictors,
        categorical_predictors = specification$categorical_predictors,
        nfolds = cfg$nfolds, seed = cfg$seed, bootstrap_B = cfg$bootstrap_B,
        thresholds = cfg$thresholds, approximate_oof_ci = cfg$approximate_oof_ci,
        parameters = cfg$parameters, h2o_memory = cfg$h2o_memory)
    }
  }
  if (!length(importance_models) && stage %in% c("all", "importance", "plot"))
    message("No supported importance model requested; Bayes has no native table. Skipping importance/plots.")
  if (length(importance_models) && stage %in% c("all", "importance")) {
    # Neither extraction function refits a model or changes thresholds.
    native <- setdiff(importance_models, "SVM")
    if (length(native)) extract_h2o_varimp(models_root, schemes = schemes, models = native)
    if ("SVM" %in% importance_models)
      extract_svm_varimp(models_root, cfg$data_file, schemes = schemes,
        evaluation_dataset = cfg$importance_dataset, n_perm = cfg$n_perm, seed = cfg$seed)
  }
  if (length(importance_models) && stage %in% c("all", "plot"))
    plot_importance_results(models_root, schemes, importance_models, evaluation_dataset = cfg$importance_dataset)
  invisible(TRUE)
}
