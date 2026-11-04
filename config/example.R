# A small, fully simulated SVM example; no Java server is needed.
source(file.path(repo_root, "config", "project.R"), local = TRUE)
config$data_file <- "data/example/model_dat.tsv"
config$labels_file <- "data/example/sample_label.tsv"
config$output_dir <- "outputs/example"
config$models <- "SVM"
config$importance_models <- "SVM"
config$schemes <- config$schemes["simplify_family_abun_model_01"]
config$nfolds <- 3L
config$bootstrap_B <- 200L
config$n_perm <- 10L
