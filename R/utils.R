# Shared path and saved-model handling; no dependency on a workstation layout.
absolute_path <- function(path, root) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) stop("Invalid path.")
  path <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) path <- file.path(root, path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

read_model_manifest <- function(development_dir) {
  manifest <- readRDS(file.path(development_dir, "manifest.rds"))
  manifest$model_files <- lapply(manifest$model_files, function(path) {
    # Prefer co-located models so copied runs never load a model from the old run.
    local <- file.path(development_dir, "models", basename(path))
    if (file.exists(local)) return(normalizePath(local, winslash = "/"))
    resolved <- absolute_path(path, development_dir)
    if (!file.exists(resolved)) stop("Saved model missing: ", local)
    resolved
  })
  manifest
}

resolve_config <- function(config_file, repo_root) {
  config_file <- absolute_path(config_file, repo_root)
  environment <- new.env(parent = .GlobalEnv)
  environment$repo_root <- normalizePath(repo_root, winslash = "/")
  sys.source(config_file, envir = environment)
  if (!exists("config", envir = environment, inherits = FALSE)) stop("Config file must define config <- list(...).")
  cfg <- environment$config
  for (name in c("data_file", "labels_file", "output_dir"))
    cfg[[name]] <- absolute_path(cfg[[name]], repo_root)
  cfg$models_root <- if (is.null(cfg$models_root)) file.path(cfg$output_dir, "model_outputs") else absolute_path(cfg$models_root, repo_root)
  cfg$importance_models <- intersect(cfg$importance_models, cfg$models)
  cfg$config_file <- config_file
  if (!length(cfg$schemes) || is.null(names(cfg$schemes)) || anyDuplicated(names(cfg$schemes)) ||
      any(!grepl("^[A-Za-z0-9_.-]+$", names(cfg$schemes))) || any(names(cfg$schemes) %in% c(".", "..")))
    stop("Schemes must have unique, simple directory names.")
  if (!length(cfg$models) || anyDuplicated(cfg$models) || !all(cfg$models %in% names(default_model_parameters()))) stop("Invalid models.")
  if (!all(cfg$importance_models %in% c("GBM", "GLM", "RF", "DL", "xGboost", "SVM")))
    stop("Importance supports GBM, GLM, RF, DL, xGboost and SVM; Bayes has no native table.")
  if (!cfg$importance_dataset %in% c("temporal", "development")) stop("Invalid importance_dataset.")
  cfg
}

read_cohorts <- function(cfg) {
  if (!file.exists(cfg$data_file) || !file.exists(cfg$labels_file))
    stop("Input files missing. Edit config/project.R or generate the simulated example first.")
  data <- read.delim(cfg$data_file, check.names = FALSE, stringsAsFactors = FALSE,
    colClasses = setNames("character", cfg$id_col))
  labels <- read.delim(cfg$labels_file, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
  label_id <- if (is.null(cfg$label_id_col)) names(labels)[1L] else cfg$label_id_col
  if (!all(c(label_id, cfg$label_col) %in% names(labels)) || !cfg$id_col %in% names(data)) stop("ID/label column missing.")
  if (anyNA(data[[cfg$id_col]]) || anyDuplicated(data[[cfg$id_col]]) || any(!nzchar(data[[cfg$id_col]])) ||
      anyNA(labels[[label_id]]) || anyDuplicated(labels[[label_id]])) stop("IDs must be nonmissing and unique.")
  if (!setequal(data[[cfg$id_col]], labels[[label_id]])) stop("Data and cohort labels must contain exactly the same patient IDs.")
  split <- labels[[cfg$label_col]][match(data[[cfg$id_col]], labels[[label_id]])]
  if (anyNA(split) || !all(split %in% c(cfg$train_label, cfg$test_label)) || identical(cfg$train_label, cfg$test_label)) stop("Invalid cohort labels.")
  train <- data[split == cfg$train_label, , drop = FALSE]
  test <- data[split == cfg$test_label, , drop = FALSE]
  for (scheme in names(cfg$schemes)) {
    specification <- cfg$schemes[[scheme]]
    if (!all(specification$predictors %in% names(data))) stop("Missing predictors for scheme: ", scheme)
    prepared <- prepare_modelFun_data(train, test, cfg$outcome_col, cfg$id_col,
      specification$predictors, specification$categorical_predictors)
    outcome <- unique(as.character(prepared$train[[cfg$outcome_col]]))
    if (length(outcome) != 2L || anyNA(outcome) || !cfg$positive %in% outcome) stop("Training outcome must contain two classes including positive.")
    negative <- setdiff(outcome, cfg$positive)
    y <- check_data(prepared$train, cfg$id_col, cfg$outcome_col, specification$predictors, cfg$positive, negative)
    check_data(prepared$test, cfg$id_col, cfg$outcome_col, specification$predictors, cfg$positive, negative)
    make_folds(y, cfg$seed, cfg$nfolds)
  }
  list(train = train, test = test)
}
