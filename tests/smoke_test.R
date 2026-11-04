#!/usr/bin/env Rscript
# Integration test with synthetic data; writes to a temporary directory only.
self <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- dirname(dirname(normalizePath(self)))
source(file.path(root, "run_pipeline.R"))
load_pipeline(root)
source(file.path(root, "scripts/make_example_data.R"))
work <- tempfile("modeling_smoke_")
dir.create(work)
make_example_data(file.path(work, "data"))
cfg <- resolve_config("config/example.R", root)
cfg$data_file <- file.path(work, "data/model_dat.tsv")
cfg$labels_file <- file.path(work, "data/sample_label.tsv")
cfg$models_root <- file.path(work, "model_outputs")
if ("--h2o" %in% commandArgs(trailingOnly = TRUE)) {
  cfg$models <- c("GBM", "SVM")
  cfg$importance_models <- cfg$models
  cfg$parameters$GBM$ntrees <- 10L
}
run_pipeline(cfg, "all")
scheme <- names(cfg$schemes)[1L]
development <- file.path(cfg$models_root, scheme, "development")
manifest <- readRDS(file.path(development, "manifest.rds"))
stopifnot(all(!grepl("^(/|[A-Za-z]:)", unlist(manifest$model_files))))
# Relocate models and recover from the moved directory, without touching the old one.
moved <- file.path(work, "relocated_development")
dir.create(moved)
stopifnot(file.copy(file.path(development, "manifest.rds"), moved))
dir.create(file.path(moved, "models"))
stopifnot(all(file.copy(list.files(file.path(development, "models"), full.names = TRUE), file.path(moved, "models"))))
loaded <- read_model_manifest(moved)
stopifnot(all(startsWith(unlist(loaded$model_files), normalizePath(moved))))
for (model in cfg$importance_models) {
  folder <- file.path(cfg$models_root, "variable_importance", scheme)
  result <- read_importance(file.path(folder, paste0(model, "_varimp.tsv")), model)
  chart <- make_importance_plot(result)
  stopifnot(nrow(chart$data) == min(10L, nrow(result)))
  stopifnot(all(file.info(file.path(folder, paste0(model, "_top10.", c("pdf", "png"))))$size > 1000))
}
# Negative and zero permutation scores must stay valid and correctly normalized.
result <- format_svm_varimp(c(a = .3, b = .1, c = -.2))
stopifnot(isTRUE(all.equal(result$percentage, c(.75, .25, 0))),
  all(format_svm_varimp(c(a = 0, b = 0))$percentage == 0),
  svm_auc(c(0L, 1L), c(.9, .1)) == 0)
cat("PASS: synthetic modeling, locked temporal evaluation, importance extraction, Top10 PDF/PNG, model relocation and score normalization.\n")
cat("Temporary test outputs:", work, "\n")
