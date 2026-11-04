#!/usr/bin/env Rscript
# Entirely simulated data. No patient records are read or embedded.
make_example_data <- function(outdir, seed = 2026L) {
  set.seed(seed)
  n_development <- 80L; n_temporal <- 28L; n <- n_development + n_temporal
  smoke <- rbinom(n, 1, .35); lung <- rbinom(n, 1, .25)
  cd8 <- rnorm(n, 25, 8); il8 <- rlnorm(n, 3, .65)
  antibody <- sample(0:2, n, replace = TRUE, prob = c(.55, .3, .15))
  virus <- rbinom(n, 1, .4) * rlnorm(n, -2, 1)
  risk <- plogis(-.5 + .6 * smoke + .8 * lung - .04 * (cd8 - 25) - .012 * (il8 - 20) + .7 * (virus > 0))
  data <- data.frame(sample_code = sprintf("SIM_%04d", seq_len(n)),
    Group = ifelse(rbinom(n, 1, risk), "CIP", "Non-CIP"),
    Smoking_history = smoke, Chronic_lung_disease_any = lung, CD8_Ts = cd8,
    IL_8 = il8, autoAb_domain_positive = antibody, CD4_Th = rnorm(n, 40, 8), IL_2R = rlnorm(n, 6, .4),
    Herpesviridae_abun = virus, Herpesviridae_det = as.integer(virus > 0), check.names = FALSE)
  for (prefix in c("Simplexvirus_humanalpha1", "Roseolovirus_humanbeta7", "Roseolovirus_humanbeta6b",
      "Roseolovirus_humanbeta6a", "Lymphocryptovirus_humangamma4", "Cytomegalovirus_humanbeta5")) {
    detected <- rbinom(n, 1, .3)
    data[[paste0(prefix, "_abun")]] <- detected * rlnorm(n, -3, 1)
    data[[paste0(prefix, "_det")]] <- detected
  }
  # Exercise training-fold imputation without changing labels or cohort identity.
  data$IL_8[c(7L, 43L, 95L)] <- NA_real_
  labels <- data.frame(sample_code = data$sample_code,
    label = c(rep("train", n_development), rep("test", n_temporal)))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  for (file in c("model_dat.tsv", "sample_label.tsv")) if (file.exists(file.path(outdir, file))) stop("Example input already exists: ", file)
  write.table(data, file.path(outdir, "model_dat.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(labels, file.path(outdir, "sample_label.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(list(data = data, labels = labels))
}
if (sys.nframe() == 0L) {
  self <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
  root <- dirname(dirname(normalizePath(self)))
  args <- commandArgs(trailingOnly = TRUE)
  target <- if (length(args)) args[1L] else file.path(root, "data", "example")
  make_example_data(target)
  cat("Simulated inputs created:", target, "\n")
}
