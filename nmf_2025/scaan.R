library(LESYMAP)
library(oro.nifti)
library(ANTsRCore)

# file paths
lesion_dir    <- 'path to directory with binary lesion masks (MNI normalized)' # should have patient IDs in filenames
output_dir    <- 'output directory path'
behavior_file <- 'path to file with IDs, AQ, demographics' # expects columns for ID, aphasia_quotient, age, gender (either 1/0 or M/F), months_post_onset, years_of_education

# behavior data
behavior_data <- read.csv(behavior_file, stringsAsFactors = FALSE)
behavior_data$id <- as.character(behavior_data$id)
behavior_data$aphasia_quotient <- as.numeric(behavior_data$aphasia_quotient)
rownames(behavior_data) <- behavior_data$id
behavior_data$gender <- ifelse(behavior_data$gender == "M", 1, 0)

behavior_vector <- behavior_data$aphasia_quotient
names(behavior_vector) <- behavior_data$id
cov_df <- behavior_data[, c("aphasia_quotient", "gender", "age", "months_post_onset", "years_of_education")]

# lesion data
lesion_files <- list.files(lesion_dir, pattern = "\\.nii$", full.names = TRUE)
lesion_ids   <- sub("\\.nii$", "", basename(lesion_files))
match_idx    <- lesion_ids %in% rownames(behavior_data)
lesion_files <- lesion_files[match_idx]
lesion_ids   <- lesion_ids[match_idx]
behavior_vector <- behavior_vector[lesion_ids]
cov_df <- cov_df[lesion_ids, , drop = FALSE]
lesion_images <- imageFileNames2ImageList(lesion_files)

# SCCAN function using lesymap
run_lesymap_sccan <- function(covariates = NULL, suffix, correctByLesSize = "none") {
  args <- list(
    lesions.list = lesion_images,
    behavior = behavior_vector,
    method = "sccan",
    correctByLesSize = correctByLesSize,
    pThreshold = 0.05,
    minSubjectPerVoxel = "10%",
    validateSparseness = TRUE,
    parallel = TRUE
  )
  
  if (!is.null(covariates)) {
    args$covariates <- covariates
    args$behavior_column <- "aphasia_quotient"
  }
  
  res <- do.call(lesymap, args)
  
  cat("SCCAN Result for", suffix, "\n")
  
  if (!is.null(res$stat.img)) {
    antsImageWrite(res$stat.img, file.path(output_dir, paste0("statMap_sccan_", suffix, ".nii.gz")))
    if (!is.null(res$rawWeights.img)) {
      antsImageWrite(res$rawWeights.img, file.path(output_dir, paste0("rawWeights_sccan_", suffix, ".nii.gz")))
    }
    cat("Optimal Sparseness: ", res$sparseness, "\n")
    cat("Predictive Correlation (cv): ", res$predictiveCor, "\n\n")
  } else {
    cat("No valid SCCAN result\n\n")
  }
  
  return(res)
}


### No covariates
result_sccan_noCov <- run_lesymap_sccan(covariates = NULL, suffix = "noCov")
antsImageWrite(result_sccan_noCov$mask.img, file.path(output_dir, "threshold_mask_sccan.nii.gz"))

### Lesion size regressed out
result_sccan_lv <- run_lesymap_sccan(covariates = NULL, suffix = "lv", correctByLesSize = "behavior")

### Lesion size + all covariates
all_covs <- cov_df[, c("gender", "age", "months_post_onset", "years_of_education")]
result_sccan_allCov <- run_lesymap_sccan(covariates = all_covs, suffix = "allCov", correctByLesSize = "behavior")

