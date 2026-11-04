# Export H2O native variable importance using the four-column reference schema.
# No training, permutation importance, or evaluation-data scoring is performed.

native_importance_map <- function(prep) {
  original <- character()
  for (nm in names(prep$specs)) {
    spec <- prep$specs[[nm]]
    count <- if (spec$type == "numeric") 1L else length(spec$levels) - 1L
    if (count) original <- c(original, rep(nm, count))
  }
  encoded <- sprintf("v%05d", seq_along(original))
  retained <- original[encoded %in% prep$columns]
  retained <- names(prep$specs)[names(prep$specs) %in% retained]
  setNames(retained, sprintf("h%05d", seq_along(retained)))
}

format_native_importance <- function(importance, prep, aggregate_categorical = FALSE) {
  columns <- c("variable", "relative_importance", "scaled_importance", "percentage")
  importance <- as.data.frame(importance, stringsAsFactors = FALSE)
  if (!all(columns %in% names(importance))) stop("Unexpected H2O importance schema.")
  importance <- importance[, columns, drop = FALSE]
  mapping <- native_importance_map(prep)
  encoded <- as.character(importance$variable)
  codes <- sub("\\..*$", "", encoded)
  original <- unname(mapping[codes])
  if (anyNA(original)) stop("Cannot map encoded variable(s): ", paste(encoded[is.na(original)], collapse = ", "))
  scores <- as.numeric(importance$relative_importance)
  if (any(!is.finite(scores)) || any(scores < 0)) stop("Invalid native importance scores.")
  if (aggregate_categorical) {
    importance <- aggregate(scores, by = list(variable = original), FUN = sum)
    names(importance)[2L] <- "relative_importance"
    maximum <- max(importance$relative_importance)
    total <- sum(importance$relative_importance)
    importance$scaled_importance <- if (maximum > 0) importance$relative_importance / maximum else 0
    importance$percentage <- if (total > 0) importance$relative_importance / total else 0
  } else {
    # Preserve the exact native scores, changing only the displayed variable name.
    labels <- original
    for (i in seq_along(encoded)) {
      suffix <- substring(encoded[i], nchar(codes[i]) + 1L)
      if (nzchar(suffix)) {
        spec <- prep$specs[[original[i]]]
        if (spec$type == "categorical" && grepl("^\\.L[0-9]+$", suffix)) {
          level <- as.integer(sub("^\\.L", "", suffix))
          if (level > length(spec$levels)) stop("Unknown category: ", encoded[i])
          labels[i] <- paste0(original[i], "=", spec$levels[level])
        } else labels[i] <- paste0(original[i], suffix)
      }
    }
    importance$variable <- labels
  }
  importance <- importance[order(-importance$relative_importance, importance$variable), columns, drop = FALSE]
  rownames(importance) <- NULL
  importance
}

extract_h2o_varimp <- function(
  models_root,
  schemes = NULL, models = c("GBM", "GLM", "RF", "DL", "xGboost"),
  outdir = file.path(models_root, "h2o_varimp"),
  aggregate_categorical = FALSE, overwrite = FALSE
) {
  if (!requireNamespace("h2o", quietly = TRUE)) stop("Install h2o first.")
  supported <- c("GBM", "GLM", "RF", "DL", "xGboost")
  if (!length(models) || anyDuplicated(models) || !all(models %in% supported)) {
    stop("Native export supports GBM, GLM, RF, DL and xGboost; Bayes/SVM have no native table in these runs.")
  }
  root <- models_root
  if (is.null(schemes)) {
    candidates <- list.dirs(root, recursive = FALSE, full.names = FALSE)
    schemes <- candidates[vapply(candidates, function(s)
      file.exists(file.path(root, s, "development/manifest.rds")), logical(1))]
  }
  if (!length(schemes) || anyDuplicated(schemes)) stop("No schemes or duplicate schemes.")
  manifests <- setNames(lapply(schemes, function(s)
    read_model_manifest(file.path(root, s, "development"))), schemes)
  for (scheme in schemes) {
    man <- manifests[[scheme]]
    if (!all(models %in% man$config$models)) stop("Model absent from scheme: ", scheme)
    if (!all(file.exists(unlist(man$model_files[models])))) stop("Model file missing: ", scheme)
  }
  if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)) && !overwrite) {
    stop("Output directory is not empty; choose another outdir or set overwrite=TRUE.")
  }
  h2o::h2o.init(max_mem_size = "4G", nthreads = 1L)
  h2o::h2o.no_progress()
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  exported <- list()
  for (scheme in schemes) {
    man <- manifests[[scheme]]
    directory <- file.path(outdir, scheme)
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    for (name in models) {
      read_native <- function() {
        existing <- h2o::h2o.ls()$key
        model <- h2o::h2o.loadModel(man$model_files[[name]])
        if (!model@model_id %in% existing) on.exit(h2o::h2o.rm(model@model_id), add = TRUE)
        # Stored native table; no replacement by a different scoring method.
        table <- h2o::h2o.varimp(model)
        if (is.null(table) || !nrow(as.data.frame(table))) stop("Native importance unavailable: ", scheme, "/", name)
        table
      }
      native <- read_native()
      importance <- format_native_importance(native, man$preprocessor, aggregate_categorical)
      file <- file.path(directory, paste0(scheme, "_", name, "_varimp.tsv"))
      write.table(importance, file, sep = "\t", row.names = FALSE, quote = FALSE,
        na = "NA", fileEncoding = "UTF-8")
      exported[[paste(scheme, name, sep = "/")]] <- importance
      message("Exported: ", scheme, "/", name, " (", nrow(importance), " rows)")
    }
  }
  readme <- c(
    "H2O native variable importance",
    "Each scheme folder contains <scheme>_<model>_varimp.tsv.",
    "Columns exactly match the reference: variable, relative_importance, scaled_importance, percentage.",
    "Rows are sorted by descending relative_importance; no row names or quoting.",
    "GBM/GLM/RF/DL/xGboost use h2o.varimp() on the saved final models; no refitting or permutation.",
    "Bayes and e1071 radial SVM have no corresponding native table here and are not exported.",
    if (aggregate_categorical)
      "Categorical levels are summed per original variable; scaled_importance and percentage are recalculated."
    else
      "Native scores and rows are preserved. Encoded categorical levels are renamed, e.g. Smoking_history=1.",
    "scaled_importance = relative_importance / maximum; percentage = relative_importance / total (0-1 fraction).",
    "Native scores have algorithm-specific definitions and are intended for ranking within a model."
  )
  writeLines(readme, file.path(outdir, "README.txt"))
  message("Complete: ", normalizePath(outdir, winslash = "/"))
  invisible(exported)
}
