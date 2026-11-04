# Paths are relative to the repository root, not the shell working directory.
# Put private input files under data/private/ (excluded by .gitignore).
source(file.path(repo_root, "config", "feature_schemes.R"), local = TRUE)
config <- list(
  data_file = "data/private/model_dat.tsv",
  labels_file = "data/private/sample_label.tsv",
  output_dir = "outputs/project",
  models_root = NULL, # Optional existing model_outputs directory, including legacy runs.
  id_col = "sample_code", outcome_col = "Group", positive = "CIP",
  label_id_col = NULL, # NULL = first column of the cohort-label table.
  label_col = "label", train_label = "train", test_label = "test",
  models = c("GBM", "GLM", "RF", "DL", "xGboost", "Bayes", "SVM"),
  importance_models = c("GBM", "GLM", "RF", "DL", "xGboost", "SVM"),
  schemes = feature_schemes(),
  nfolds = 5L, seed = 2026L, bootstrap_B = 2000L,
  thresholds = NULL, # Select Youden threshold on development OOF predictions.
  approximate_oof_ci = FALSE, h2o_memory = "4G",
  parameters = default_model_parameters(),
  importance_dataset = "temporal", n_perm = 100L
)
