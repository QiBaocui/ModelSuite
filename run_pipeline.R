#!/usr/bin/env Rscript
# Usage: Rscript run_pipeline.R [config/project.R] [all|check|model|importance|plot]
load_pipeline <- function(repo_root, envir = .GlobalEnv) {
  for (file in c("utils.R", "model_5fold_temporal.R", "importance_h2o.R",
      "importance_svm.R", "plot_importance.R", "pipeline.R"))
    source(file.path(repo_root, "R", file), local = envir, encoding = "UTF-8")
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  argument <- grep("^--file=", commandArgs(), value = TRUE)
  repo_root <- dirname(normalizePath(sub("^--file=", "", argument[1L]), winslash = "/"))
  args <- commandArgs(trailingOnly = TRUE)
  if (any(args %in% c("--help", "-h"))) {
    cat("Usage: Rscript run_pipeline.R [config/project.R] [all|check|model|importance|plot]\n")
    quit(status = 0L)
  }
  if (getRversion() < "4.1.0") stop("R >= 4.1 is required.")
  load_pipeline(repo_root)
  config_file <- if (length(args)) args[1L] else "config/project.R"
  stage <- if (length(args) >= 2L) args[2L] else "all"
  cfg <- resolve_config(config_file, repo_root)
  run_pipeline(cfg, stage)
}
