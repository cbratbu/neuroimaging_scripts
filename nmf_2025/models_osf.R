# File paths
h_data_path <- "path/to/H_data.csv"       # H_data file with patient loadings from NMF output; column names should be 'Atom_X'. No ID column is expected, so row order should match ID order in the other tables.
tissue_data_path <- "path/to/tissue_data.csv"  # Table of lesion burden data per ROI; must contain an "ID" column.
covariate_data_path <- "path/to/covariate_data.csv"  # Contains ID column, aphasia_quotient, lesion_volume, gender, age, months_post_onset, years_of_education

library(dplyr)

# Load data files
H_data <- read.csv(h_data_path, stringsAsFactors = FALSE)
tissue_data <- read.csv(tissue_data_path, stringsAsFactors = FALSE)
covariate_data <- read.csv(covariate_data_path, stringsAsFactors = FALSE)

# Create single data table
data_demo <- cbind(covariate_data, H_data) 
data_full <- merge(data_demo, tissue_data, by = "ID", sort = FALSE)

atom_columns <- colnames(H_data)
k <- length(atom_columns)

### Single-Predictor Regressions on all variables
predictor_vars <- setdiff(colnames(data_full), c("aphasia_quotient", "ID"))
single_var_models <- list()
for (var in predictor_vars) {
  formula_str <- paste("aphasia_quotient ~", var)
  single_var_models[[var]] <- lm(as.formula(formula_str), data = data_full)
}

### Full Models 

# function for LOOCV
calc_loocv_rmse <- function(model_formula, data) {
  n <- nrow(data)
  predictions <- numeric(n)
  for (i in 1:n) {
    train <- data[-i, ]
    test <- data[i, ]
    mod <- lm(as.formula(model_formula), data = train)
    predictions[i] <- predict(mod, newdata = test)
  }
  actual <- data$aphasia_quotient
  sqrt(mean((actual - predictions)^2, na.rm = TRUE))
}

# Baseline model (lesion volume and demographics)
baseline_formula <- "aphasia_quotient ~ gender + age + months_post_onset + years_of_education + lesion_volume"
baseline_model <- lm(as.formula(baseline_formula), data = data_full)
baseline_model$loocv_rmse <- calc_loocv_rmse(baseline_formula, data_full)

# single atom models
single_atom_models <- list()
for (atom in atom_columns) {
  formula_str <- paste("aphasia_quotient ~ gender + age + months_post_onset + years_of_education + lesion_volume *", atom)
  model_obj <- lm(as.formula(formula_str), data = data_full)
  model_obj$loocv_rmse <- calc_loocv_rmse(formula_str, data_full)
  single_atom_models[[atom]] <- model_obj
}

# model with all atoms
all_atoms <- paste(atom_columns, collapse = " + ")
full_formula <- paste0("aphasia_quotient ~ gender + age + months_post_onset + years_of_education + lesion_volume * (", all_atoms, ")")
full_model <- lm(as.formula(full_formula), data = data_full)
full_model$loocv_rmse <- calc_loocv_rmse(full_formula, data_full)
