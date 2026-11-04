# Organize existing importance tables and plot independent Top10 figures.
# Horizontal Top10 charts with the established blue gradient style.
# Does not train models or recompute importance.

read_importance <- function(path, model) {
  df <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  expected <- c("variable", "relative_importance", "scaled_importance", "percentage")
  if (!identical(names(df), expected) || !nrow(df) || anyNA(df) ||
      any(!nzchar(df$variable)) || anyDuplicated(df$variable)) stop("Invalid table: ", path)
  for (column in expected[-1L]) {
    if (!is.numeric(df[[column]]) || any(!is.finite(df[[column]])))
      stop("Invalid numeric scores: ", path)
  }
  tolerance <- 1e-7
  score <- df$relative_importance
  maximum <- if (model == "SVM") max(abs(score)) else max(score)
  expected_scaled <- if (maximum > 0) score / maximum else rep(0, nrow(df))
  positive <- pmax(score, 0)
  expected_percentage <- if (sum(positive) > 0) positive / sum(positive) else rep(0, nrow(df))
  if (model != "SVM" && any(score < 0)) stop("Unexpected negative native scores: ", path)
  if (max(abs(df$scaled_importance - expected_scaled)) > tolerance ||
      max(abs(df$percentage - expected_percentage)) > tolerance)
    stop("Score normalization mismatch: ", path)
  df
}

make_importance_plot <- function(df) {
  df <- df[order(-df$scaled_importance, seq_len(nrow(df))), , drop = FALSE]
  df <- head(df, 10L)
  df$Feature <- factor(df$variable, levels = rev(df$variable))
  df$Fill <- ifelse(df$scaled_importance < 0, NA_real_, df$scaled_importance)
  df$Label_hjust <- ifelse(df$scaled_importance < 0, 1.18, -0.18)
  lower <- min(c(0, df$scaled_importance))
  lower <- if (lower < 0) lower - .16 else 0
  base_theme <- ggplot2::theme_classic(base_size = 12, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15, margin = ggplot2::margin(b = 7)),
      plot.subtitle = ggplot2::element_text(color = "grey40", size = 10, margin = ggplot2::margin(b = 12)),
      axis.text = ggplot2::element_text(color = "grey15"),
      axis.title = ggplot2::element_text(face = "bold", size = 12),
      plot.margin = ggplot2::margin(14, 18, 14, 14))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = scaled_importance, y = Feature)) +
    ggplot2::geom_col(ggplot2::aes(fill = Fill), width = .68, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", scaled_importance), hjust = Label_hjust),
      size = 3.3, color = "grey20") +
    ggplot2::scale_fill_gradient(low = "#AED6F1", high = "#154360", limits = c(0, 1), na.value = "#C87545") +
    ggplot2::scale_x_continuous(limits = c(lower, 1.15),
      breaks = round(seq(if (lower < 0) ceiling(lower * 5) / 5 else 0, 1, .2), 1),
      expand = ggplot2::expansion(mult = 0)) +
    ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 32),
      expand = ggplot2::expansion(add = .6)) +
    ggplot2::labs(title = NULL, subtitle = NULL, x = "Scaled Importance", y = NULL) +
    base_theme + ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.major.x = ggplot2::element_line(color = "grey90", linetype = "dashed", linewidth = .3),
      axis.ticks.y = ggplot2::element_blank())
  if (any(df$scaled_importance < 0))
    p <- p + ggplot2::geom_vline(xintercept = 0, color = "grey45", linewidth = .4)
  list(plot = p, data = df)
}

plot_importance_results <- function(models_root, schemes, models,
  evaluation_dataset = "temporal", outdir = file.path(models_root, "variable_importance"),
  overwrite = FALSE) {
  for (package in c("ggplot2", "stringr"))
    if (!requireNamespace(package, quietly = TRUE)) stop("Install package: ", package)
  if (utils::packageVersion("ggplot2") < "3.4.0") stop("ggplot2 >= 3.4.0 is required.")
  if (!capabilities("cairo")) stop("Cairo support is required for PDF and PNG.")
  if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)) && !overwrite)
    stop("Plot output exists; use a fresh output directory.")
  tables <- list(); sources <- list(); predictors <- list()
  for (scheme in schemes) {
    man <- read_model_manifest(file.path(models_root, scheme, "development"))
    predictors[[scheme]] <- man$config$predictors
    for (model in models) {
      source_root <- if (model == "SVM") file.path(models_root, "svm_varimp", evaluation_dataset) else file.path(models_root, "h2o_varimp")
      path <- file.path(source_root, scheme, paste0(scheme, "_", model, "_varimp.tsv"))
      key <- paste(scheme, model, sep = "/")
      tables[[key]] <- read_importance(path, model)
      sources[[key]] <- path
    }
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  for (scheme in schemes) {
    directory <- file.path(outdir, scheme)
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    for (model in models) {
      key <- paste(scheme, model, sep = "/")
      target <- file.path(directory, paste0(model, "_varimp.tsv"))
      if (!file.copy(sources[[key]], target, overwrite = overwrite)) stop("Cannot copy table: ", key)
      if (unname(tools::md5sum(sources[[key]])) != unname(tools::md5sum(target))) stop("Copy mismatch.")
      figure <- make_importance_plot(tables[[key]])
      prefix <- file.path(directory, paste0(model, "_top10"))
      ggplot2::ggsave(paste0(prefix, ".pdf"), figure$plot, device = grDevices::cairo_pdf,
        width = 9.2, height = 5.8, units = "in", bg = "white")
      ggplot2::ggsave(paste0(prefix, ".png"), figure$plot, device = "png", type = "cairo",
        width = 9.2, height = 5.8, units = "in", dpi = 300, bg = "white")
      message("Plotted: ", key, "; Top", nrow(figure$data))
    }
  }
  explanation <- c("Variable importance tables and independent Top10 PDF/PNG charts.",
    paste("Models:", paste(models, collapse = ", ")),
    "Tables retain all rows; charts display at most 10 rows sorted by scaled_importance (descending).",
    "Native H2O scores preserve expanded categorical levels. SVM permutes original categorical variables jointly.",
    paste0("SVM permutation evaluation dataset: ", evaluation_dataset, ". This is not OOF importance."),
    "Negative SVM scores are retained and plotted in orange. Other bars use a blue gradient.",
    "percentage is a 0-1 fraction. Native and permutation scores are not directly comparable across algorithms.",
    "Bayes has no native importance table and is excluded from this output.", "", "Scheme predictor lists:")
  for (scheme in schemes) explanation <- c(explanation, paste0(scheme, ": ", paste(predictors[[scheme]], collapse = ", ")))
  writeLines(explanation, file.path(outdir, "README.txt"))
  invisible(tables)
}
