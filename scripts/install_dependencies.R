#!/usr/bin/env Rscript
# Only run this installer deliberately; pipeline scripts never install packages.
if (getRversion() < "4.1.0") stop("R >= 4.1 is required.")
packages <- c("h2o", "e1071", "pROC", "PRROC", "ggplot2", "stringr")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
if (utils::packageVersion("ggplot2") < "3.4.0") stop("Upgrade ggplot2 to >= 3.4.0.")
cat("Dependencies available. For saved H2O models, match the original H2O client/server version.\n")
