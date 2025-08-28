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

# VLSM function using lesymap
run_lesymap <- function(covariates = NULL, suffix, correctByLesSize = "none", 
                        patchinfo = NA, multipleComparison = "FWERperm") {
  args <- list(
    lesions.list = lesion_images,
    behavior = behavior_vector,
    method = "BMfast",
    correctByLesSize = correctByLesSize,
    multipleComparison = multipleComparison,
    pThreshold = 0.05,
    minSubjectPerVoxel = "10%",
    nperm = 1000,
    patchinfo = patchinfo,
    parallel = TRUE
  )
  
  if (!is.null(covariates)) {
    args$covariates <- covariates
    args$behavior_column <- "aphasia_quotient"
  }
  
  res <- do.call(lesymap, args)
  antsImageWrite(res$stat.img, file.path(output_dir, paste0("statMap_", suffix, ".nii.gz")))
  
  if (!is.null(res$perm.FWEthresh)) {
    cat("FWE threshold for", suffix, ":", res$perm.FWEthresh, "\n")
  }
  
  return(res)
}

### No covariates

# Uncorrected
result_noCov_uncorr <- run_lesymap(covariates = NULL,
                                   suffix = "uncorr_noCov",
                                   correctByLesSize = "none",
                                   multipleComparison = "none")

# Save shared outputs
patchData <- result_noCov_uncorr$patchinfo
antsImageWrite(result_noCov_uncorr$mask.img, file.path(output_dir, "threshold_mask.nii.gz"))
antsImageWrite(result_noCov_uncorr$average.img, file.path(output_dir, "average_lesion_map.nii.gz"))

# FWER
result_noCov_fwer <- run_lesymap(covariates = NULL,
                                 suffix = "fwer_noCov",
                                 correctByLesSize = "none",
                                 patchinfo = patchData,
                                 multipleComparison = "FWERperm")

# FDR
result_noCov_fdr <- run_lesymap(covariates = NULL,
                                suffix = "fdr_noCov",
                                correctByLesSize = "none",
                                patchinfo = patchData,
                                multipleComparison = "FDR")

# Bonferroni
result_noCov_bonf <- run_lesymap(covariates = NULL,
                                 suffix = "bonf_noCov",
                                 correctByLesSize = "none",
                                 patchinfo = patchData,
                                 multipleComparison = "Bonferroni")


### Lesion size regressed out

# Uncorrected
result_lv_uncorr <- run_lesymap(covariates = NULL,
                                suffix = "uncorr_lv",
                                correctByLesSize = "behavior",
                                patchinfo = patchData,
                                multipleComparison = "none")

# FWER
result_lv_fwer <- run_lesymap(covariates = NULL,
                              suffix = "fwer_lv",
                              correctByLesSize = "behavior",
                              patchinfo = patchData,
                              multipleComparison = "FWERperm")

# FDR
result_lv_fdr <- run_lesymap(covariates = NULL,
                             suffix = "fdr_lv",
                             correctByLesSize = "behavior",
                             patchinfo = patchData,
                             multipleComparison = "FDR")

# Bonferroni
result_lv_bonf <- run_lesymap(covariates = NULL,
                              suffix = "bonf_lv",
                              correctByLesSize = "behavior",
                              patchinfo = patchData,
                              multipleComparison = "Bonferroni")


### Lesion size + all covariates

all_covs <- cov_df[, c("gender", "age", "months_post_onset", "years_of_education")]

# Uncorrected
result_allCov_uncorr <- run_lesymap(covariates = all_covs,
                                    suffix = "uncorr_allCov",
                                    correctByLesSize = "behavior",
                                    patchinfo = patchData,
                                    multipleComparison = "none")

# FWER
result_allCov_fwer <- run_lesymap(covariates = all_covs,
                                  suffix = "fwer_allCov",
                                  correctByLesSize = "behavior",
                                  patchinfo = patchData,
                                  multipleComparison = "FWERperm")

# FDR
result_allCov_fdr <- run_lesymap(covariates = all_covs,
                                 suffix = "fdr_allCov",
                                 correctByLesSize = "behavior",
                                 patchinfo = patchData,
                                 multipleComparison = "FDR")

# Bonferroni
result_allCov_bonf <- run_lesymap(covariates = all_covs,
                                  suffix = "bonf_allCov",
                                  correctByLesSize = "behavior",
                                  patchinfo = patchData,
                                  multipleComparison = "Bonferroni")

