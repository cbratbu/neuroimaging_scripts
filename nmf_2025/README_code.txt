Please contact the Center for Brain Recovery at Boston University (brainrec@bu.edu) or the corresponding author, Emerson Kropp (ekropp@bu.edu), with any questions or concerns related to these data.

Code used for preparation of spared tissue data, NMF decomposition, and linear regression analysis:


gray_matter_spared.py -> Python code to calculate percent of spared tissue in given ROIs (used with AAL3 atlas). Inputs are directory containing binary lesion mask .nii files, directory containing ROI .nii files, and output directory. Output is a .csv file listing percent spared tissue per ROI (column) per subject (row). **The tractotron feature in BCBtoolbox (Foulon et al., 2018) was used for calculating lesion overlap with white matter tracts; this data was concatenated with spared gray matter data to build the complete spared tissue data file used in analysis. 

kim_tidor_k_selection.py -> Python code comparing results of dimensionality reduction under different methods (NMF, SVD, SVD on random matrix) to determine k value to use in NMF (based on method from Kin & Tidor, 2003). Takes input .csv containing region-wise tissue data per subject, outputs figure plotting accuracy of each method and reconstruction errors per k value. 

nmf_script.py -> Python code to run NMF decomposition of region-wise tissue data. Inputs are .csv file of region-wise tissue data per subject as well as pre-selected k value, outputs are NMF matrices W and H. 

nmf_aq_models.R -> R code to run various linear regression models comparing results using NMF atom patient loadings with various other predictors. Expects NMF matrix H results, as well as table containing WAB AQ, age, gender (either M/F or 1/0), months_post_onset, years_of_education, an lesion_volume, with corresponding patient IDs. Can also input region-wise tissue data to run models on each individual ROI predicting AQ. Runs linear regressions predicting AQ using only lesion volume + demographics, LV/demographics with each NMF atom individually, and LV/demographics with all NMF atoms included as predictors. 

vlsm_lesymap.R -> R code to run univariate VLSM using LESYMAP package (Pustina et al., 2018). Inputs are binarized lesion maps, WAB AQ data, and demographic covariates. Runs VLSM with no covariates, lesion volume regressed out, or LV + all demographics regressed out. Each is output with no correction, FDR correction, Bonferroni FWER correction, or permutation-based FWER correction. Also outputs map of lesion overlap per subject and threshold map of voxels lesioned in >10% of subjects. 

sccan.R -> R code to run multivariate SCCAN using LESYMAP package (Pustina et al., 2018). Inputs are binarized lesion maps, WAB AQ data, and demographic covariates. Runs SCCAN with no covariates, lesion volume regressed out, or LV + all demographics regressed out. Uses internal cross-validation to determine optimal sparseness, if a value is identified it outputs a statistical map.

