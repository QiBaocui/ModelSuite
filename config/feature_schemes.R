# Exact predictor and categorical-variable lists for the 12 original schemes.
feature_schemes <- function() {
  categorical_clinical <- c("Smoking_history", "Chronic_lung_disease_any", "autoAb_domain_positive")
  species <- c("Simplexvirus_humanalpha1", "Roseolovirus_humanbeta7", "Roseolovirus_humanbeta6b",
    "Roseolovirus_humanbeta6a", "Lymphocryptovirus_humangamma4", "Cytomegalovirus_humanbeta5")
  schemes <- list()
  for (clinical_type in c("simplify_01", "simplify_02", "overall_01")) {
    clinical <- if (clinical_type == "overall_01")
      c("Smoking_history", "Chronic_lung_disease_any", "CD8_Ts", "IL_8", "CD4_Th", "IL_2R", "autoAb_domain_positive") else
      c("Smoking_history", "Chronic_lung_disease_any", "CD8_Ts",
        if (clinical_type == "simplify_01") "IL_8" else "IL_2R", "autoAb_domain_positive")
    for (level in c("family", "species")) for (measurement in c("abun", "det")) {
      virus <- paste0(if (level == "family") "Herpesviridae" else species, "_", measurement)
      name <- paste(if (clinical_type == "overall_01") "overall" else "simplify", level, measurement,
        "model", if (clinical_type == "simplify_02") "02" else "01", sep = "_")
      schemes[[name]] <- list(predictors = c(clinical, virus),
        categorical_predictors = c(categorical_clinical, if (measurement == "det") virus else character()))
    }
  }
  schemes
}
