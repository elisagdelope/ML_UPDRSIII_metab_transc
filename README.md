# Table of contents
* [Introduction](#introduction)
* [Content](#content)
* [Data](#data)
* [Requirements](#requirements)
* [License](#license)

## Introduction
This repository contains the code for ML analyses performed in Chapter 5 of my PhD thesis "Interpretable machine learning on omics data for the study of UPDRS III prognosis". The project consists on predicting the Unified Parkinson’s Disease Rating Scale Part III (UPDRS III) motor scores (mild/severe when classification) from whole blood transcriptomics and blood plasma metabolomics using measurements from the baseline clinical visit, and temporal or dynamic features engineered from a short temporal series of 4 and 3 timepoints, respectively, from the PPMI cohort and the LuxPARK cohort, aiming at identifying molecular and higher-level functional fingerprints linked specifically to the motor symptoms in PD.

## Content
The code covers the following main tasks and analyses:

1. Preprocessing of the datasets, including unsupervised filters.
2. Classification of mild vs. severe UPDRS III using baseline measurements in a nested cross-validation setting with LASSO feature selection and evaluation on test set, for molecular and higher-level functional representations as predictors.
3. Classification of mild vs. severe UPDRS III using temporal (dynamic) features in a nested cross-validation setting with LASSO feature selection, for molecular and higher-level functional representations as predictors.
4. Regression of UPDRS III scores, both using baseline and temporal features in a nested cross-validation setting with LASSO feature selection, for molecular and higher-level functional representations as predictors.
5. Feature importance with SHAP values.
6. Post-hoc analysis with Friedman tests to compare performance results, and plots on the performance.

There is a README.md inside each directory with corresponding explanation for each script:

*ppmi_analyses* contains code related to analysis on transcriptomics data from PPMI.

*luxpark_analyses* contains code related to analysis on metabolomics data from LuxPARK.


## Data
The public transcriptomics data used in this project was derived from the Parkinson’s Progression Markers Initiative (https://www.ppmi-info.org/, RNAseq - IR3).
The metabolomics data from LuxPARK is not publicly available as it is linked to the Luxembourg Parkinson’s Study and its internal regulations. Any requests for accessing the dataset can be directed to request.ncer-pd@uni.lu.

## Requirements
The code for ML modeling was implemented in Python (3.9.13), and that for the pre-processing steps, was implemented in R (R 4.0.3). It has been tested on both current Mac (Ventura) and Linux operating systems (Rocky Linux 8.7 (Green Obsidian)), relying on multiple Python, R, and BioConductor packages correspondingly that are listed at the beginning of each script. The code should be compatible with later versions of Python and R installed on current Mac, Linux or Windows systems. **For the ML modeling, it is necessary to download and install the NestedCV() class implemented in the digipd_ml package and environment from repository [digipd_ml](https://gitlab.lcsb.uni.lu/elisa.gomezdelope/digipd_ml).**

## License
TBD


