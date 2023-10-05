# ma-cc-rare
Code for The Impact of Corrections Methods on Rare-Event Meta-Analysis

This folder contains code to reproduce the simulation study (Section 4 in the paper) and the case study example (Section 5 in the paper). 

## For the simulation study:
Create two directories: RDSDataSets and Results. Make sure the following packages are installed: meta, metafor, BiasedUrn, future.apply, and nleqslv.

<ins>0 Miscellaneous scripts needed for subsequent scripts (do not run):</ins>
- 0_DGMs.R
  - Contains binary event meta-analysis data generating functions
- 0_CCCode.R
  - Contains code to implement various continuity corrections to individual MA data sets
  - Calls 0_MAResults.R
- 0_ApplyAllCCs.R
  - Applies various continuity corrections to individual MA data sets
  - Calls 0_CCCode.R
- 0_MAResults.R
  - Runs an individual MA data set through metabin/metafor, separated by heterogeneity variance estimator
  - Calls 0_IPM.R
- 0_IPM.R
  - Contains code to implement the improved Paule-Mandel (IPM) heterogeneity variance estimator

<ins>1 Scripts to generate meta-analysis data:</ins>
- 1_GenerateMAData.R
  - Generates binary event meta-analysis data
  - Data generated is stored in the folder RDSDataSets
  - Need to set the appropriate amount of workers 
  - Calls 0_DGMs.R
- 1_GetpCFixedResults.R
  - Applies various continuity corrections, heterogeneity estimators, and pooling methods to previously generated data sets
  - This code is for pCFixed – do the same things for pRandomB
  - Need to set the appropriate amount of workers 
  - Calls 0_ApplyAllCCs.R and 0_MAResults.R
- 1_GetpCFixedGLMMResults.R
  - Applies the GLMM to previously generated data sets
  - This code is for pCFixed – do the same things for pRandomB
  - Need to set the appropriate amount of workers 
Results can then be summarized in terms of Type I error rate, relative power, confidence interval coverage, median confidence interval width, mean squared error, and median bias, as outlined in the paper.

## For the case study example:
Make sure the following packages are installed: tidyverse, readxl, and meta.
1. Anti-TNF.xlsx
  - Anti-TNF data set
2. Case_Study_Anti-TNF.R 
  - Applies various continuity corrections, heterogeneity estimators, and pooling methods to the Anti-TNF data set
  - Reads in Anti-TNF.xlsx
  - Calls 0_ApplyAllCCs.R and 0_MAResults.R (defined above)
