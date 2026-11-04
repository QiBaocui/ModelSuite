# One shared stratified K-fold CV (default 5) + locked temporal validation.
# Requires R >= 4.1, h2o, e1071, pROC, PRROC and Java for the H2O server.
# No automatic hyperparameter tuning, ensemble, or model selection.
# Default thresholds maximise Youden on development OOF predictions only.
# See README.md and config/example.R before running.

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = TRUE, row.names = FALSE,
              na = "NA", fileEncoding = "UTF-8")
}

check_packages <- function() {
  pkgs <- c("h2o", "e1071", "pROC", "PRROC")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
}

check_data <- function(d, id_col, outcome_col, predictors, positive, negative) {
  if (!is.data.frame(d) || anyDuplicated(names(d))) stop("Input must have unique column names.")
  required <- c(id_col, outcome_col, predictors)
  if (!all(required %in% names(d))) stop("Missing columns: ", paste(setdiff(required, names(d)), collapse = ", "))
  if (!length(predictors) || anyDuplicated(predictors) ||
      any(predictors %in% c(id_col, outcome_col))) stop("Invalid predictor list.")
  if (anyNA(d[[id_col]]) || anyDuplicated(as.character(d[[id_col]]))) stop("Patient IDs must be nonmissing and unique.")
  labels <- as.character(d[[outcome_col]])
  if (anyNA(labels) || !all(labels %in% c(negative, positive))) stop("Outcome labels do not match positive/negative.")
  y <- as.integer(labels == positive)
  if (length(unique(y)) != 2L) stop("Both outcome classes are required.")
  ok <- vapply(d[predictors], function(z) is.numeric(z) || is.factor(z) || is.character(z) || is.logical(z), logical(1))
  if (!all(ok)) stop("Unsupported predictor types: ", paste(predictors[!ok], collapse = ", "))
  y
}

check_nfolds <- function(nfolds) {
  if (!is.numeric(nfolds) || length(nfolds) != 1L || !is.finite(nfolds) ||
      nfolds < 2 || nfolds > .Machine$integer.max || nfolds != floor(nfolds)) {
    stop("nfolds must be a single integer >= 2.")
  }
  as.integer(nfolds)
}

make_folds <- function(y, seed, nfolds = 5L) {
  nfolds <- check_nfolds(nfolds)
  if (min(table(y)) < nfolds) stop("Each outcome class needs at least ", nfolds,
    " patients for stratified ", nfolds, "-fold CV.")
  set.seed(seed)
  fold <- integer(length(y))
  for (g in 0:1) {
    ii <- which(y == g)
    fold[ii] <- sample(rep(seq_len(nfolds), length.out = length(ii)))
  }
  fold
}

# All preprocessing is fitted ONLY on the current training subset.
# Numeric: training median imputation; all-missing numeric predictors removed.
# Categorical: training mode imputation; K-1 dummy coding; unseen levels ->
# training mode WITH a warning. This rule is fixed before validation.
# Constant encoded columns removed; remaining columns centred/scaled.
fit_preprocessor <- function(d, y, feature_selector = NULL) {
  selected <- names(d)
  if (!is.null(feature_selector)) selected <- feature_selector(d, y)
  if (!is.character(selected) || !length(selected) || anyDuplicated(selected) ||
      !all(selected %in% names(d))) stop("feature_selector must return valid unique predictor names.")
  specs <- list()
  for (nm in selected) {
    z <- d[[nm]]
    if (is.numeric(z)) {
      z[!is.finite(z)] <- NA_real_
      if (all(is.na(z))) next
      specs[[nm]] <- list(type = "numeric", fill = median(z, na.rm = TRUE))
    } else {
      z <- as.character(z)
      if (all(is.na(z))) next
      counts <- sort(table(z), decreasing = TRUE)
      mode <- names(counts)[1L]
      levels <- c(mode, sort(setdiff(unique(z[!is.na(z)]), mode)))
      specs[[nm]] <- list(type = "categorical", fill = mode, levels = levels)
    }
  }
  if (!length(specs)) stop("No usable predictors in this training subset.")
  prep <- list(specs = specs)
  raw <- encode_data(prep, d, warn_unseen = FALSE)
  keep <- vapply(raw, function(z) length(unique(z)) > 1L, logical(1))
  if (!any(keep)) stop("All encoded predictors are constant in a training fold.")
  prep$columns <- names(raw)[keep]
  prep$centre <- vapply(raw[keep], mean, numeric(1))
  prep$scale <- vapply(raw[keep], sd, numeric(1))
  prep
}

encode_data <- function(prep, d, warn_unseen = TRUE) {
  columns <- list()
  for (nm in names(prep$specs)) {
    spec <- prep$specs[[nm]]
    z <- d[[nm]]
    if (spec$type == "numeric") {
      if (!is.numeric(z)) stop("Numeric predictor changed type: ", nm)
      z[!is.finite(z)] <- spec$fill
      columns[[length(columns) + 1L]] <- z
    } else {
      z <- as.character(z)
      unseen <- !is.na(z) & !z %in% spec$levels
      if (warn_unseen && any(unseen)) warning("Unseen levels in ", nm, ": ",
        paste(unique(z[unseen]), collapse = ", "), "; mapped to training mode: ", spec$fill)
      z[is.na(z) | unseen] <- spec$fill
      for (lv in spec$levels[-1L]) columns[[length(columns) + 1L]] <- as.numeric(z == lv)
    }
  }
  if (!length(columns)) return(data.frame(row.names = seq_len(nrow(d))))
  names(columns) <- sprintf("v%05d", seq_along(columns))
  as.data.frame(columns, check.names = FALSE)
}

apply_preprocessor <- function(prep, d) {
  raw <- encode_data(prep, d)[prep$columns]
  raw[] <- Map(function(z, a, b) (z - a) / b, raw, prep$centre, prep$scale)
  raw
}

# H2O receives native categorical predictors, in particular so Naive Bayes
# uses categorical likelihoods rather than treating dummy indicators as Gaussian.
# SVM receives the numeric dummy matrix from apply_preprocessor().
apply_native_preprocessor <- function(prep, d) {
  mapping <- preprocessor_map(prep)
  columns <- list()
  for (nm in names(prep$specs)) {
    s <- prep$specs[[nm]]
    retained <- mapping[mapping$original == nm & mapping$retained, , drop = FALSE]
    if (!nrow(retained)) next
    z <- d[[nm]]
    if (s$type == "numeric") {
      if (!is.numeric(z)) stop("Numeric predictor changed type: ", nm)
      z[!is.finite(z)] <- s$fill
      code <- retained$encoded[1L]
      z <- (z - prep$centre[[code]]) / prep$scale[[code]]
    } else {
      z <- as.character(z)
      unseen <- !is.na(z) & !z %in% s$levels
      if (any(unseen)) warning("Unseen levels in ", nm, "; mapped to training mode: ", s$fill)
      z[is.na(z) | unseen] <- s$fill
      z <- factor(match(z, s$levels), levels = seq_along(s$levels),
                  labels = paste0("L", seq_along(s$levels)))
    }
    columns[[nm]] <- z
  }
  # Safe, stable column names within this fitted preprocessor.
  names(columns) <- sprintf("h%05d", seq_along(columns))
  as.data.frame(columns, check.names = FALSE)
}

preprocessor_map <- function(prep) {
  rows <- list()
  for (nm in names(prep$specs)) {
    s <- prep$specs[[nm]]
    if (s$type == "numeric") {
      rows[[length(rows) + 1L]] <- data.frame(original = nm, encoding = "numeric", fill = as.character(s$fill))
    } else for (lv in s$levels[-1L]) {
      rows[[length(rows) + 1L]] <- data.frame(original = nm,
        encoding = paste0("level=", lv, "; reference=", s$levels[1L]), fill = s$fill)
    }
  }
  ans <- do.call(rbind, rows)
  ans$encoded <- sprintf("v%05d", seq_len(nrow(ans)))
  ans$retained <- ans$encoded %in% prep$columns
  ans
}

# Fixed illustrative settings, NOT claimed to be optimal for any cohort size.
# Change before the analysis; tuning requires a separate nested-CV design.
default_model_parameters <- function() list(
  GBM = list(ntrees = 100L, max_depth = 2L, min_rows = 5,
             learn_rate = 0.03, stopping_rounds = 0L),
  GLM = list(family = "binomial", alpha = 0, lambda = 0.01,
             lambda_search = FALSE, standardize = FALSE),
  RF = list(ntrees = 500L, max_depth = 5L, min_rows = 3),
  DL = list(hidden = c(8L), epochs = 50, activation = "Rectifier",
            l2 = 0.01, standardize = FALSE, stopping_rounds = 0L,
            reproducible = TRUE, train_samples_per_iteration = -1L),
  xGboost = list(ntrees = 100L, max_depth = 2L, min_rows = 3,
                 learn_rate = 0.03, stopping_rounds = 0L, backend = "cpu"),
  Bayes = list(laplace = 1, min_sdev = 0.01),
  SVM = list(kernel = "radial", cost = 1)
)

fit_engine <- function(name, X, y, seed, parameters) {
  set.seed(seed)
  if (name == "SVM") {
    args <- c(list(x = as.matrix(X), y = factor(y, levels = 0:1),
                   scale = FALSE, probability = TRUE, type = "C-classification"), parameters)
    return(do.call(e1071::svm, args))
  }
  frame <- h2o::as.h2o(data.frame(X, .outcome = factor(y, levels = 0:1)))
  on.exit(h2o::h2o.rm(frame), add = TRUE)
  fun <- switch(name, GBM = h2o::h2o.gbm, GLM = h2o::h2o.glm,
    RF = h2o::h2o.randomForest, DL = h2o::h2o.deeplearning,
    xGboost = h2o::h2o.xgboost, Bayes = h2o::h2o.naiveBayes)
  # Manual shared outer folds; nfolds=0 avoids a second CV inside each fit.
  do.call(fun, c(list(x = names(X), y = ".outcome", training_frame = frame,
                     nfolds = 0L, seed = seed), parameters))
}

predict_engine <- function(name, model, X) {
  if (name == "SVM") {
    pred <- predict(model, newdata = as.matrix(X), probability = TRUE)
    probs <- attr(pred, "probabilities")
    if (is.null(probs) || !"1" %in% colnames(probs)) stop("Missing SVM positive probability column.")
    p <- probs[, "1"]
  } else {
    frame <- h2o::as.h2o(X)
    on.exit(h2o::h2o.rm(frame), add = TRUE)
    pred_frame <- h2o::h2o.predict(model, frame)
    on.exit(h2o::h2o.rm(pred_frame), add = TRUE)
    probs <- as.data.frame(pred_frame)
    if (!"p1" %in% names(probs)) stop("Missing H2O p1 column; expected factor levels 0,1.")
    p <- probs$p1
  }
  if (length(p) != nrow(X) || any(!is.finite(p)) || any(p < 0 | p > 1)) stop("Invalid predicted probabilities.")
  as.numeric(p)
}

safe_glm <- function(formula, data) {
  bad <- FALSE
  fit <- tryCatch(withCallingHandlers(glm(formula, data = data, family = binomial()),
    warning = function(w) { bad <<- TRUE; invokeRestart("muffleWarning") }), error = function(e) NULL)
  if (is.null(fit) || bad || !isTRUE(fit$converged) || any(!is.finite(coef(fit)))) return(NULL)
  fit
}

calibration_statistics <- function(y, p) {
  lp <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  ans <- c(Calibration_intercept_CITL = NA_real_,
           Calibration_slope = NA_real_, Calibration_joint_intercept = NA_real_)
  if (length(unique(y)) != 2L) return(ans)
  # CITL is intercept with slope fixed at 1, not the joint-fit intercept.
  a <- safe_glm(y ~ 1 + offset(lp), data.frame(y, lp))
  if (!is.null(a)) ans[1L] <- unname(coef(a)[1L])
  separated <- max(lp[y == 0]) <= min(lp[y == 1]) ||
    max(lp[y == 1]) <= min(lp[y == 0])
  if (sd(lp) > 1e-8 && !separated) {
    b <- safe_glm(y ~ lp, data.frame(y, lp))
    if (!is.null(b)) ans[2:3] <- unname(coef(b)[c(2L, 1L)])
  }
  ans
}

basic_statistics <- function(y, p) {
  if (length(unique(y)) != 2L) return(c(AUC = NA_real_, PR_AUC = NA_real_,
    Brier = mean((p - y)^2), calibration_statistics(y, p)))
  roc <- pROC::roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)
  pr <- PRROC::pr.curve(scores.class0 = p[y == 1], scores.class1 = p[y == 0], curve = FALSE)
  c(AUC = as.numeric(pROC::auc(roc)), PR_AUC = pr$auc.integral,
    Brier = mean((p - y)^2), calibration_statistics(y, p))
}

smooth_calibration <- function(y, p, grid) {
  if (length(unique(y)) != 2L || length(unique(p)) < 5L) return(rep(NA_real_, length(grid)))
  # Low-complexity flexible logistic spline; descriptive in tiny samples.
  lp <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  fit <- safe_glm(y ~ splines::ns(lp, df = 2), data.frame(y, lp))
  if (is.null(fit)) return(rep(NA_real_, length(grid)))
  tryCatch(as.numeric(predict(fit, newdata = data.frame(
    lp = qlogis(pmin(pmax(grid, 1e-6), 1 - 1e-6))), type = "response")),
    error = function(e) rep(NA_real_, length(grid)))
}

net_benefit <- function(y, p, thresholds) {
  vapply(thresholds, function(t) mean(y == 1 & p >= t) -
    mean(y == 0 & p >= t) * t / (1 - t), numeric(1))
}

ci_quantiles <- function(z, min_valid) {
  z <- z[is.finite(z)]
  if (length(z) < min_valid) return(c(NA_real_, NA_real_))
  as.numeric(quantile(z, c(0.025, 0.975), names = FALSE))
}

# Candidates are finite probabilities in [0,1], consistent with p >= threshold.
# Include endpoints; ties in Youden (within 1e-12) use the lowest threshold.
select_oof_youden <- function(y, p) {
  if (length(y) != length(p) || anyNA(y) || !all(y %in% 0:1) ||
      length(unique(y)) != 2L || any(!is.finite(p)) || any(p < 0 | p > 1)) {
    stop("Youden selection requires valid binary outcomes and probabilities.")
  }
  candidates <- sort(unique(c(0, p, 1)))
  bins <- match(p, candidates)
  tp <- rev(cumsum(rev(tabulate(bins[y == 1], nbins = length(candidates)))))
  fp <- rev(cumsum(rev(tabulate(bins[y == 0], nbins = length(candidates)))))
  sensitivity <- tp / sum(y == 1)
  specificity <- 1 - fp / sum(y == 0)
  youden <- sensitivity + specificity - 1
  best <- which(youden >= max(youden) - 1e-12)
  j <- best[1L]
  data.frame(threshold = candidates[j], Youden = youden[j],
    OOF_sensitivity = sensitivity[j], OOF_specificity = specificity[j],
    tied_candidates = length(best), tie_rule = "lowest_threshold_tolerance_1e-12",
    source = "development_oof_youden", classification_rule = "probability >= threshold")
}

check_threshold_settings <- function(thresholds, threshold_method) {
  threshold_method <- match.arg(threshold_method, c("oof_youden", "fixed"))
  if (threshold_method == "oof_youden") {
    if (!is.null(thresholds)) stop("For oof_youden, leave thresholds=NULL.")
  } else if (!is.numeric(thresholds) || !length(thresholds) ||
      any(!is.finite(thresholds)) || any(thresholds <= 0 | thresholds >= 1)) {
    stop("For fixed thresholds, provide probability cutoff(s) in (0,1).")
  }
  threshold_method
}

# Read per-model thresholds, with compatibility for older fixed-threshold runs.
locked_thresholds <- function(cfg, model_name) {
  if (is.null(cfg$thresholds_by_model)) {
    if (identical(cfg$threshold_method, "oof_youden")) stop("Missing locked OOF thresholds.")
    return(cfg$thresholds)
  }
  ans <- cfg$thresholds_by_model[[model_name]]
  if (!length(ans) || any(!is.finite(ans)) || any(ans < 0 | ans > 1)) {
    stop("Invalid locked threshold for ", model_name)
  }
  ans
}

threshold_table <- function(y, p, thresholds, uncertainty) {
  ratio_ci <- function(a, n) {
    if (!n) return(c(NA_real_, NA_real_, NA_real_))
    ci <- if (uncertainty) binom.test(a, n)$conf.int else c(NA_real_, NA_real_)
    c(a / n, ci)
  }
  do.call(rbind, lapply(thresholds, function(t) {
    positive <- p >= t
    TP <- sum(positive & y == 1); FP <- sum(positive & y == 0)
    FN <- sum(!positive & y == 1); TN <- sum(!positive & y == 0)
    se <- ratio_ci(TP, TP + FN); sp <- ratio_ci(TN, TN + FP)
    ppv <- ratio_ci(TP, TP + FP); npv <- ratio_ci(TN, TN + FN)
    data.frame(threshold = t, TP, FP, FN, TN,
      Sensitivity = se[1], Sensitivity_low = se[2], Sensitivity_high = se[3],
      Specificity = sp[1], Specificity_low = sp[2], Specificity_high = sp[3],
      PPV = ppv[1], PPV_low = ppv[2], PPV_high = ppv[3],
      NPV = npv[1], NPV_low = npv[2], NPV_high = npv[3],
      Accuracy = (TP + TN) / length(y),
      F1 = if ((2 * TP + FP + FN) > 0) 2 * TP / (2 * TP + FP + FN) else NA_real_)
  }))
}

draw_plots <- function(y, p, stats, cal, dca, path, title) {
  # One figure per file; path is the output prefix (normally "evaluation").
  plots <- list(
    roc = function() {
      roc <- pROC::roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)
      plot(1 - roc$specificities, roc$sensitivities, type = "l", col = "#2369a3",
        xlab = "False positive rate", ylab = "Sensitivity", xlim = c(0, 1), ylim = c(0, 1),
        main = sprintf("ROC: AUC = %.3f", stats["AUC"]))
      abline(0, 1, lty = 2, col = "grey")
    },
    pr = function() {
      pr <- PRROC::pr.curve(scores.class0 = p[y == 1], scores.class1 = p[y == 0], curve = TRUE)
      plot(pr$curve[, 1], pr$curve[, 2], type = "l", xlim = c(0, 1), ylim = c(0, 1),
        col = "#2369a3", xlab = "Recall", ylab = "Precision",
        main = sprintf("PR: AUC = %.3f", stats["PR_AUC"]))
      abline(h = mean(y), lty = 2, col = "grey")
    },
    calibration = function() {
      plot(c(0, 1), c(0, 1), type = "n", xlab = "Predicted probability",
        ylab = "Observed probability", main = "Calibration (flexible spline)")
      if (all(is.finite(cal$low)) && all(is.finite(cal$high))) polygon(
        c(cal$predicted, rev(cal$predicted)), c(cal$low, rev(cal$high)),
        col = adjustcolor("#2369a3", alpha.f = 0.18), border = NA)
      abline(0, 1, lty = 2, col = "grey")
      lines(cal$predicted, cal$observed, col = "#2369a3", lwd = 2)
      rug(p[y == 0], side = 1, col = "#777777")
      rug(p[y == 1], side = 3, col = "#2369a3")
      if (!any(is.finite(cal$observed))) text(.5, .5, "Curve not estimable")
    },
    decision = function() {
      upper <- max(c(dca$net_benefit, dca$treat_all, 0), na.rm = TRUE)
      lower <- max(-0.15, min(c(dca$net_benefit, dca$treat_all, 0), na.rm = TRUE))
      plot(dca$threshold, dca$net_benefit, type = "l", col = "#2369a3", lwd = 2,
        ylim = c(lower - 0.02, upper + 0.02), xlab = "Threshold probability",
        ylab = "Net benefit", main = "Decision curve (exploratory)")
      lines(dca$threshold, dca$treat_all, lty = 2, col = "#b64e31")
      abline(h = 0, lty = 3)
      legend("topright", c("Model", "Treat all", "Treat none"),
        col = c("#2369a3", "#b64e31", "black"), lty = 1:3, bty = "n", cex = .8)
    }
  )
  draw <- function(draw_plot) {
    old <- par(mfrow = c(1, 1), mar = c(4.5, 4.5, 4, 1))
    on.exit(par(old))
    draw_plot()
    mtext(title, side = 3, line = 2.5, cex = .9)
  }
  for (name in names(plots)) {
    plot_path <- paste0(path, "_", name)
    pdf(paste0(plot_path, ".pdf"), width = 6, height = 6)
    tryCatch(draw(plots[[name]]), finally = dev.off())
    png(paste0(plot_path, ".png"), width = 1800, height = 1800, res = 300)
    tryCatch(draw(plots[[name]]), finally = dev.off())
  }
}

# CI bootstraps below resample evaluation patients/predictions; NEVER refit models.
# They do not constitute repeated CV or bootstrap internal validation.
# OOF CI disabled by default: overlapping training folds invalidate naive inference.
evaluate_predictions <- function(y, p, ids, model_name, dataset, outdir,
  thresholds, dca_thresholds, bootstrap_B = 2000L, seed = 1L,
  approximate_oof_ci = FALSE, fold = rep(NA_integer_, length(y)),
  threshold_source = "prespecified_fixed") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  uncertainty <- dataset == "temporal" || (dataset == "internal_oof" && approximate_oof_ci)
  scope <- if (dataset == "temporal") "conditional_on_locked_model" else if (uncertainty)
    "approximate_fixed_OOF_predictions_NOT_full_pipeline_CI" else "not_estimated"
  s <- basic_statistics(y, p)
  metrics <- data.frame(model = model_name, dataset = dataset, N = length(y),
    events = sum(y), prevalence = mean(y), metric = names(s), estimate = as.numeric(s),
    CI_low = NA_real_, CI_high = NA_real_, CI_method = "not_estimated",
    CI_scope = scope, bootstrap_valid = NA_integer_)
  grid <- seq(max(min(p), 1e-6), min(max(p), 1 - 1e-6), length.out = 101L)
  cal <- data.frame(predicted = grid, observed = smooth_calibration(y, p, grid),
    low = NA_real_, high = NA_real_, bootstrap_valid = NA_integer_)
  dca <- data.frame(threshold = dca_thresholds,
    net_benefit = net_benefit(y, p, dca_thresholds),
    treat_all = mean(y) - (1 - mean(y)) * dca_thresholds / (1 - dca_thresholds),
    treat_none = 0, low = NA_real_, high = NA_real_, bootstrap_valid = NA_integer_)
  if (uncertainty) {
    set.seed(seed)
    bs <- matrix(NA_real_, bootstrap_B, length(s), dimnames = list(NULL, names(s)))
    bc <- matrix(NA_real_, bootstrap_B, nrow(cal))
    bd <- matrix(NA_real_, bootstrap_B, nrow(dca))
    for (b in seq_len(bootstrap_B)) {
      ii <- sample.int(length(y), replace = TRUE)
      bs[b, ] <- basic_statistics(y[ii], p[ii])
      bc[b, ] <- smooth_calibration(y[ii], p[ii], grid)
      bd[b, ] <- net_benefit(y[ii], p[ii], dca_thresholds)
    }
    min_valid <- ceiling(0.8 * bootstrap_B)
    for (j in seq_along(s)) {
      ci <- ci_quantiles(bs[, j], min_valid)
      metrics$CI_low[j] <- ci[1]; metrics$CI_high[j] <- ci[2]
      metrics$bootstrap_valid[j] <- sum(is.finite(bs[, j]))
      metrics$CI_method[j] <- "patient_pairs_bootstrap_percentile"
    }
    if (dataset == "temporal") {
      roc <- pROC::roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)
      ci <- tryCatch(as.numeric(pROC::ci.auc(roc, method = "delong")), error = function(e) NULL)
      if (!is.null(ci)) {
        metrics$CI_low[1] <- ci[1]; metrics$CI_high[1] <- ci[3]
        metrics$CI_method[1] <- "DeLong_conditional_on_locked_model"
        metrics$bootstrap_valid[1] <- NA_integer_
      }
    }
    ci_cal <- apply(bc, 2, ci_quantiles, min_valid = min_valid)
    cal$low <- ci_cal[1, ]; cal$high <- ci_cal[2, ]
    cal$bootstrap_valid <- colSums(is.finite(bc))
    ci_dca <- apply(bd, 2, ci_quantiles, min_valid = min_valid)
    dca$low <- ci_dca[1, ]; dca$high <- ci_dca[2, ]
    dca$bootstrap_valid <- colSums(is.finite(bd))
    if (any(colSums(is.finite(bs)) < min_valid)) warning(model_name, " ", dataset,
      ": some CIs unavailable: <80% valid bootstrap estimates; inspect metrics.")
  }
  write_tsv(data.frame(ID = ids, y = y, probability = p, fold = fold), file.path(outdir, "predictions.tsv"))
  write_tsv(metrics, file.path(outdir, "metrics.tsv"))
  selected_on_evaluation <- dataset == "internal_oof" && threshold_source == "development_oof_youden"
  threshold_uncertainty <- uncertainty && !selected_on_evaluation
  th <- threshold_table(y, p, thresholds, threshold_uncertainty)
  th$threshold_source <- threshold_source
  th$interpretation <- if (selected_on_evaluation) "exploratory_threshold_selected_on_same_OOF_data" else
    if (dataset == "apparent") "apparent_training_performance" else "evaluation_at_locked_threshold"
  th$CI_scope <- if (threshold_uncertainty) scope else "not_estimated"
  th$CI_method <- if (threshold_uncertainty) "Clopper_Pearson_exact_conditional" else "not_estimated"
  write_tsv(th, file.path(outdir, "threshold_metrics.tsv"))
  write_tsv(cal, file.path(outdir, "calibration_curve.tsv"))
  write_tsv(dca, file.path(outdir, "decision_curve.tsv"))
  roc <- pROC::roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)
  write_tsv(data.frame(FPR = 1 - roc$specificities, TPR = roc$sensitivities), file.path(outdir, "roc_curve.tsv"))
  pr <- PRROC::pr.curve(scores.class0 = p[y == 1], scores.class1 = p[y == 0], curve = TRUE)
  write_tsv(data.frame(recall = pr$curve[, 1], precision = pr$curve[, 2]), file.path(outdir, "pr_curve.tsv"))
  draw_plots(y, p, s, cal, dca, file.path(outdir, "evaluation"), paste(model_name, dataset))
  metrics
}

fit_development <- function(development, id_col, outcome_col, predictors,
  positive, negative, thresholds = NULL, dca_thresholds, outdir,
  models = c("GBM", "GLM", "RF", "DL", "xGboost", "Bayes", "SVM"),
  parameters = default_model_parameters(), seed = 2026L,
  feature_selector = NULL, approximate_oof_ci = FALSE, bootstrap_B = 2000L,
  expected_n = NULL, h2o_memory = "4G", nfolds = 5L,
  threshold_method = if (is.null(thresholds)) "oof_youden" else "fixed") {
  nfolds <- check_nfolds(nfolds)
  threshold_method <- check_threshold_settings(thresholds, threshold_method)
  check_packages()
  if (missing(positive) || missing(negative) || identical(positive, negative) ||
      length(positive) != 1L || length(negative) != 1L || anyNA(c(positive, negative))) stop("Explicit positive and negative labels required.")
  if (missing(dca_thresholds) || !is.numeric(dca_thresholds) || length(dca_thresholds) < 2L ||
      any(!is.finite(dca_thresholds)) || any(dca_thresholds <= 0 | dca_thresholds >= 1) ||
      anyDuplicated(dca_thresholds) || is.unsorted(dca_thresholds, strictly = TRUE)) stop("Prespecify increasing DCA thresholds in (0,1).")
  if (!length(models) || anyDuplicated(models) || !all(models %in% names(default_model_parameters()))) stop("Invalid models.")
  if (!all(models %in% names(parameters))) stop("Missing model parameters.")
  reserved <- c("x", "y", "training_frame", "validation_frame", "nfolds", "fold_column",
    "fold_assignment", "seed", "model_id", "probability", "scale", "type")
  if (any(vapply(parameters[models], function(z) any(names(z) %in% reserved), logical(1)))) stop("Reserved model parameter overridden.")
  if (bootstrap_B < 200L) stop("Use at least 200 CI bootstrap replicates; 2000 recommended.")
  if (!is.null(expected_n) && nrow(development) != expected_n) stop("Unexpected development cohort size.")
  y <- check_data(development, id_col, outcome_col, predictors, positive, negative)
  ids <- as.character(development[[id_col]])
  fold <- make_folds(y, seed, nfolds)
  # Use a fresh directory so old results cannot be mistaken for a complete run.
  if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE))) stop("outdir must be new or empty.")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  outdir <- normalizePath(outdir, winslash = "/", mustWork = TRUE)
  if (any(models != "SVM")) h2o::h2o.init(max_mem_size = h2o_memory, nthreads = 1L)
  if ("xGboost" %in% models && !h2o::h2o.xgboost.available()) stop(
    "H2O XGBoost backend unavailable. Use a supported Linux H2O server; model is NOT silently skipped.")
  write_tsv(data.frame(ID = ids, y = y, fold = fold), file.path(outdir, "fold_assignment.tsv"))
  cfg <- list(id_col = id_col, outcome_col = outcome_col, predictors = predictors,
    positive = positive, negative = negative, thresholds = thresholds, threshold_method = threshold_method,
    dca_thresholds = dca_thresholds, models = models, parameters = parameters,
    seed = seed, nfolds = nfolds, bootstrap_B = bootstrap_B, approximate_oof_ci = approximate_oof_ci)
  saveRDS(cfg, file.path(outdir, "analysis_config.rds"))
  oof <- matrix(NA_real_, nrow(development), length(models), dimnames = list(NULL, models))
  fold_metrics <- list()
  for (f in seq_len(nfolds)) {
    message("Outer fold ", f, "/", nfolds)
    tr <- which(fold != f); va <- which(fold == f)
    prep <- fit_preprocessor(development[tr, predictors, drop = FALSE], y[tr], feature_selector)
    saveRDS(prep, file.path(outdir, paste0("preprocessor_fold", f, ".rds")))
    write_tsv(preprocessor_map(prep), file.path(outdir, paste0("encoding_fold", f, ".tsv")))
    Xtr <- if ("SVM" %in% models) apply_preprocessor(prep, development[tr, predictors, drop = FALSE]) else NULL
    Xva <- if ("SVM" %in% models) apply_preprocessor(prep, development[va, predictors, drop = FALSE]) else NULL
    Htr <- if (any(models != "SVM")) apply_native_preprocessor(prep, development[tr, predictors, drop = FALSE]) else NULL
    Hva <- if (any(models != "SVM")) apply_native_preprocessor(prep, development[va, predictors, drop = FALSE]) else NULL
    for (m in models) {
      message("  ", m)
      mod <- fit_engine(m, if (m == "SVM") Xtr else Htr, y[tr], seed + f, parameters[[m]])
      p <- predict_engine(m, mod, if (m == "SVM") Xva else Hva)
      oof[va, m] <- p
      s <- basic_statistics(y[va], p)
      fold_metrics[[length(fold_metrics) + 1L]] <- data.frame(model = m, fold = f,
        N = length(va), events = sum(y[va]), AUC = s["AUC"], PR_AUC = s["PR_AUC"], Brier = s["Brier"])
      if (m != "SVM") h2o::h2o.rm(mod@model_id)
    }
  }
  if (anyNA(oof)) stop("Incomplete OOF predictions; no performance summary produced.")
  # Select once on development OOF predictions; never select on temporal data.
  threshold_source <- if (threshold_method == "oof_youden") "development_oof_youden" else "prespecified_fixed"
  threshold_rows <- lapply(models, function(m) {
    row <- if (threshold_method == "oof_youden") select_oof_youden(y, oof[, m]) else
      data.frame(threshold = thresholds, Youden = NA_real_, OOF_sensitivity = NA_real_,
        OOF_specificity = NA_real_, tied_candidates = NA_integer_, tie_rule = "not_applicable",
        source = threshold_source, classification_rule = "probability >= threshold")
    cbind(model = m, row)
  })
  threshold_selection <- do.call(rbind, threshold_rows)
  cfg$thresholds_by_model <- setNames(lapply(models, function(m)
    threshold_selection$threshold[threshold_selection$model == m]), models)
  cfg$threshold_selection <- threshold_selection
  write_tsv(threshold_selection, file.path(outdir, "selected_thresholds.tsv"))
  saveRDS(cfg, file.path(outdir, "analysis_config.rds"))
  write_tsv(do.call(rbind, fold_metrics), file.path(outdir, "fold_metrics_descriptive.tsv"))
  summary <- list()
  for (m in models) summary[[m]] <- evaluate_predictions(y, oof[, m], ids, m, "internal_oof",
    file.path(outdir, "internal_oof", m), locked_thresholds(cfg, m), dca_thresholds, bootstrap_B,
    seed, approximate_oof_ci, fold, threshold_source = threshold_source)
  write_tsv(do.call(rbind, summary), file.path(outdir, "internal_oof_summary.tsv"))
  # Fit one final model per algorithm using all development patients.
  prep <- fit_preprocessor(development[predictors], y, feature_selector)
  saveRDS(prep, file.path(outdir, "preprocessor_final.rds"))
  write_tsv(preprocessor_map(prep), file.path(outdir, "encoding_final.tsv"))
  X <- if ("SVM" %in% models) apply_preprocessor(prep, development[predictors]) else NULL
  H <- if (any(models != "SVM")) apply_native_preprocessor(prep, development[predictors]) else NULL
  dir.create(file.path(outdir, "models"))
  paths <- list(); apparent <- list()
  for (m in models) {
    message("Final full-development fit: ", m)
    features <- if (m == "SVM") X else H
    mod <- fit_engine(m, features, y, seed, parameters[[m]])
    p <- predict_engine(m, mod, features)
    apparent[[m]] <- evaluate_predictions(y, p, ids, m, "apparent",
      file.path(outdir, "apparent", m), locked_thresholds(cfg, m), dca_thresholds, bootstrap_B, seed,
      threshold_source = threshold_source)
    if (m == "SVM") {
      paths[[m]] <- file.path(outdir, "models", "SVM.rds")
      saveRDS(mod, paths[[m]])
    } else {
      paths[[m]] <- h2o::h2o.saveModel(mod, path = file.path(outdir, "models"), force = TRUE)
      h2o::h2o.rm(mod@model_id)
    }
  }
  write_tsv(do.call(rbind, apparent), file.path(outdir, "apparent_summary.tsv"))
  manifest <- list(config = cfg, model_files = lapply(paths, function(p) file.path("models", basename(p))), preprocessor = prep,
    development_ids = ids, fold = fold)
  saveRDS(manifest, file.path(outdir, "manifest.rds"))
  capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo_development.txt"))
  writeLines("Development complete. Lock primary model BEFORE inspecting temporal outcomes.",
             file.path(outdir, "DEVELOPMENT_COMPLETE.txt"))
  invisible(manifest)
}

# Separate stage: commit primary model BEFORE opening/evaluating temporal data.
lock_primary_model <- function(development_dir, primary_model, rationale) {
  man <- readRDS(file.path(development_dir, "manifest.rds"))
  if (!length(primary_model) || anyDuplicated(primary_model) ||
      !all(primary_model %in% man$config$models)) stop("Invalid prespecified model set.")
  if (missing(rationale) || !nzchar(rationale)) stop("Record selection rationale using development data only.")
  lock <- file.path(development_dir, "primary_model_lock.rds")
  if (file.exists(lock)) stop("Primary model already locked; will not overwrite.")
  info <- list(primary_model = primary_model, rationale = rationale, locked_at = as.character(Sys.time()))
  saveRDS(info, lock)
  write_tsv(as.data.frame(info), file.path(development_dir, "primary_model_lock.tsv"))
  invisible(info)
}

validate_temporal <- function(temporal, development_dir, outdir,
  evaluate_all_models = TRUE, expected_n = NULL, h2o_memory = "4G") {
  check_packages()
  lockfile <- file.path(development_dir, "primary_model_lock.rds")
  if (!file.exists(lockfile)) stop("Lock the primary model before temporal validation.")
  lock <- readRDS(lockfile)
  man <- read_model_manifest(development_dir); cfg <- man$config
  y <- check_data(temporal, cfg$id_col, cfg$outcome_col, cfg$predictors, cfg$positive, cfg$negative)
  if (!is.null(expected_n) && nrow(temporal) != expected_n) stop("Unexpected temporal cohort size.")
  ids <- as.character(temporal[[cfg$id_col]])
  if (any(ids %in% man$development_ids)) stop("Patient ID overlap between development and temporal cohorts.")
  if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE))) stop("Use a new or empty temporal outdir.")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  models <- if (evaluate_all_models) cfg$models else lock$primary_model
  if (any(models != "SVM")) h2o::h2o.init(max_mem_size = h2o_memory, nthreads = 1L)
  X <- if ("SVM" %in% models) apply_preprocessor(man$preprocessor, temporal[cfg$predictors]) else NULL
  H <- if (any(models != "SVM")) apply_native_preprocessor(man$preprocessor, temporal[cfg$predictors]) else NULL
  result <- list()
  for (m in models) {
    message("Temporal validation: ", m)
    mod <- if (m == "SVM") readRDS(man$model_files[[m]]) else h2o::h2o.loadModel(man$model_files[[m]])
    p <- predict_engine(m, mod, if (m == "SVM") X else H)
    result[[m]] <- evaluate_predictions(y, p, ids, m, "temporal", file.path(outdir, m),
      locked_thresholds(cfg, m), cfg$dca_thresholds, cfg$bootstrap_B, cfg$seed,
      threshold_source = if (identical(cfg$threshold_method, "oof_youden")) "development_oof_youden" else "prespecified_fixed")
    result[[m]]$analysis_role <- if (m %in% lock$primary_model) "prespecified_primary" else "secondary_exploratory"
    if (m != "SVM") h2o::h2o.rm(mod@model_id)
  }
  summary <- do.call(rbind, result)
  write_tsv(summary, file.path(outdir, "temporal_summary.tsv"))
  capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo_temporal.txt"))
  writeLines("Temporal validation complete. No refitting/recalibration/threshold optimisation performed.",
             file.path(outdir, "TEMPORAL_COMPLETE.txt"))
  invisible(summary)
}

# Convenient single-call interface for already prepared R data.frames.
# Positive class is explicit. By default thresholds are selected on development
# OOF predictions; independent validation outcomes never determine thresholds.
prepare_modelFun_data <- function(train.R, test.R, outcome_col, id_col = NULL,
  predictors = NULL, categorical_predictors = character(0)) {
  if (!is.data.frame(train.R) || !is.data.frame(test.R)) stop("train.R and test.R must be data.frames.")
  if (anyDuplicated(names(train.R)) || anyDuplicated(names(test.R))) stop("Duplicate column names.")
  if (!outcome_col %in% names(train.R) || !outcome_col %in% names(test.R)) stop("Outcome column missing.")
  if (is.null(id_col)) {
    candidates <- intersect(c("ID", "sample_code"), names(train.R))
    if (length(candidates) > 1L) stop("Both ID and sample_code present: specify id_col explicitly.")
    if (length(candidates)) id_col <- candidates[1L]
  }
  if (is.null(id_col)) {
    id_col <- ".modelFun_patient_id"
    if (id_col %in% union(names(train.R), names(test.R))) stop("Reserved generated ID column already present.")
    automatic <- function(d) identical(rownames(d), as.character(seq_len(nrow(d))))
    if (automatic(train.R) || automatic(test.R)) {
      # Positional row names cannot establish patient identity across cohorts.
      # Generate only technical row IDs, clearly recording the limitation.
      warning("No explicit patient ID and at least one cohort has automatic row names. ",
              "Generating cohort-prefixed row IDs; patient independence CANNOT be verified. ",
              "Provide id_col for real patient-overlap checking.")
      train.R[[id_col]] <- sprintf("DEV_%04d", seq_len(nrow(train.R)))
      test.R[[id_col]] <- sprintf("TIME_%04d", seq_len(nrow(test.R)))
    } else {
      train.R[[id_col]] <- rownames(train.R)
      test.R[[id_col]] <- rownames(test.R)
    }
  }
  if (!id_col %in% names(train.R) || !id_col %in% names(test.R)) stop("ID column missing in one cohort.")
  if (id_col == outcome_col) stop("ID cannot be the outcome column.")
  if (is.null(predictors)) predictors <- setdiff(names(train.R), c(id_col, outcome_col))
  if (!length(predictors) || !all(predictors %in% names(test.R))) stop("Empty predictors or predictor missing in test.R.")
  if (!all(categorical_predictors %in% predictors)) stop("categorical_predictors must be predictor names.")
  for (nm in categorical_predictors) {
    train.R[[nm]] <- factor(train.R[[nm]])
    test.R[[nm]] <- factor(test.R[[nm]])
  }
  mismatch <- vapply(predictors, function(nm)
    is.numeric(train.R[[nm]]) != is.numeric(test.R[[nm]]), logical(1))
  if (any(mismatch)) stop("Numeric/categorical type mismatch between cohorts: ",
                         paste(predictors[mismatch], collapse = ", "))
  list(train = train.R, test = test.R, id_col = id_col, predictors = predictors)
}

modelFun <- function(train.R, test.R, outdir, feature = "model",
  feature_type_list = c("GBM", "GLM", "RF", "DL", "xGboost", "Bayes", "SVM"),
  nfolds = 5L, positive = NULL, thresholds = NULL,
  dca_thresholds = seq(0.01, 0.99, by = 0.01),
  outcome_col = "Group", id_col = NULL, predictors = NULL,
  categorical_predictors = character(0),
  parameters = default_model_parameters(), seed = 2026L,
  bootstrap_B = 2000L, approximate_oof_ci = FALSE,
  feature_selector = NULL, h2o_memory = "4G",
  threshold_method = if (is.null(thresholds)) "oof_youden" else "fixed") {
  if (!is.data.frame(train.R) || !is.data.frame(test.R)) stop("train.R and test.R must be data.frames.")
  nfolds <- check_nfolds(nfolds)
  if (length(feature) != 1L || is.na(feature) || !grepl("^[[:alnum:]_.-]+$", feature) ||
      feature %in% c(".", "..")) stop("feature must be a simple directory name.")
  if (length(positive) != 1L || anyNA(positive)) stop("Specify positive='your actual positive label'.")
  threshold_method <- check_threshold_settings(thresholds, threshold_method)
  # Cohort sizes are unrestricted; class and patient-ID checks still apply.
  data <- prepare_modelFun_data(train.R, test.R, outcome_col, id_col,
    predictors, categorical_predictors)
  positive <- as.character(positive)
  labels <- unique(as.character(data$train[[outcome_col]]))
  if (anyNA(labels) || length(labels) != 2L || !positive %in% labels) stop("Training outcome must contain two classes, including positive.")
  negative <- setdiff(labels, positive)
  y <- check_data(data$train, data$id_col, outcome_col, data$predictors, positive, negative)
  check_data(data$test, data$id_col, outcome_col, data$predictors, positive, negative)
  if (min(table(y)) < nfolds) stop("Each development outcome class must have at least ", nfolds, " patients.")
  if (any(as.character(data$train[[data$id_col]]) %in% as.character(data$test[[data$id_col]]))) stop("Patient ID overlap across cohorts.")
  root <- file.path(outdir, feature)
  # Avoid expensive development if existing temporal output would prevent completion.
  for (part in c("development", "temporal")) {
    dest <- file.path(root, part)
    if (dir.exists(dest) && length(list.files(dest, all.files = TRUE, no.. = TRUE))) stop("Existing output: ", dest, "; use a new outdir or feature.")
  }
  message("Positive class = ", positive, "; negative class = ", negative,
          "; threshold method = ", threshold_method)
  message("One shared stratified ", nfolds, "-fold CV on ", nrow(data$train), " patients; final models use all ", nrow(data$train), "; temporal validation uses ", nrow(data$test), ".")
  message("DCA grid is a plotting domain; interpret only clinically meaningful thresholds.")
  development_dir <- file.path(root, "development")
  temporal_dir <- file.path(root, "temporal")
  manifest <- fit_development(data$train, data$id_col, outcome_col, data$predictors,
    positive, negative, thresholds, dca_thresholds, development_dir,
    models = feature_type_list, parameters = parameters, seed = seed,
    feature_selector = feature_selector, approximate_oof_ci = approximate_oof_ci,
    bootstrap_B = bootstrap_B, h2o_memory = h2o_memory, nfolds = nfolds,
    threshold_method = threshold_method)
  # Record ALL requested algorithms before evaluation, without selecting a winner
  # or requiring a separate manual lock/approval step.
  lock_primary_model(development_dir, feature_type_list,
    "All algorithms prespecified in modelFun() before temporal evaluation; no test-based selection.")
  temporal_summary <- validate_temporal(data$test, development_dir, temporal_dir,
    evaluate_all_models = TRUE, h2o_memory = h2o_memory)
  result <- list(
    internal_validation = read.delim(file.path(development_dir, "internal_oof_summary.tsv")),
    apparent = read.delim(file.path(development_dir, "apparent_summary.tsv")),
    temporal_validation = temporal_summary,
    selected_thresholds = manifest$config$threshold_selection,
    fold_assignment = read.delim(file.path(development_dir, "fold_assignment.tsv")),
    manifest = manifest, output_dir = normalizePath(root, winslash = "/"))
  message("Completed. Results: ", result$output_dir)
  invisible(result)
}
