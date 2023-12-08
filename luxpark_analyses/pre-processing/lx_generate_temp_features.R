# Title: lx_generate_temp_features.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script generates temporal features from gene expression data
# Usage (metabolite level): Rscript lx_generate_temp_features.R 
# Usage (aggregated level): Rscript lx_generate_temp_features.R -l pw -s mean
# Data: data from expression counts at several timepoints (>=3), temporal features functions.


# Packages ---------------------------------------------------------------------
library(readr)
library(dplyr)
library(vroom)
library(tidyr)
library(tibble)
library(stringr)
library(argparser, quietly = TRUE)
library(matrixStats)



# I/O --------------------------------------------------------------------------
analysis_name <- "02-pred-TS-UPDRS3"
IN_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
IN_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level")
OUT_DIR <- paste0(IN_DIR, "/temporal_data")
OUT_DIR_PATHWAY <- paste0(IN_DIR_PATHWAY, "/temporal_data")


# Add command line arguments
p <- arg_parser("Generate temporal features from metabolite abundance data", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "--level", help = "level of features (pw / metab)", default = "metab", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--stat", help = "aggregated statistic (mean / median / sd / pathifier / pca)", default = "mean", type = "string", nargs = '+')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
e_level = toupper(argv$level) # gene (g) / aggregations
st = tolower(argv$stat) # stat at aggregation level

if (e_level == "METAB") {
  METAB.FILE <- file.path(IN_DIR, "log_transformed_data_fte.tsv")
  features_varname <- "METABOLITES_ID"
} else { # aggregations
  METAB.FILE <- file.path(IN_DIR_PATHWAY, paste0("log_transformed_PW_", st, ".tsv"))
  features_varname <- "PATHWAY_NAME"
  
  if ((!e_level %in% c("PW", "METAB")) | (!st %in% c("mean", "median", "sd", "pathifier", "pca")) | (!file.exists(METAB.FILE))) { 
    stop("Adequate arguments were not provided")
  }
}



# Main -------------------------------------------------------------------------
if ((!dir.exists(OUT_DIR)) | (!dir.exists(OUT_DIR_PATHWAY))) {
  dir.create(OUT_DIR, recursive = T)
  dir.create(OUT_DIR_PATHWAY, recursive = T)
}



# Data reformatting ------------------------------------------------------------

# remove patients that DON'T have data for 3 timepoints (V0, V1, V2)
metab <- vroom(METAB.FILE, col_types = cols()) 
n_tp = 3
patients_list <- metab %>%
  group_by(PATIENT_ID) %>% 
  filter(n() == n_tp) %>%   
  ungroup() %>% 
  dplyr::select(PATIENT_ID) %>%
  distinct() %>%
  pull()

metab <- metab %>% 
  filter(PATIENT_ID %in% patients_list)

# calculate dynamic features from metaboomics data
source("func_temporal_features.R")
for (temp_feature in c("sd", "lm-time")) {
  if (temp_feature == "sd") {
    metab_temp <- exprs_sd(metab)
  } else if (temp_feature == "lm-time") {
    metab_temp <- exprs_lm(metab, features_varname)
  } 
  
  # remove those features having NA for the dynamic features -----
  metab_temp <- metab_temp %>%
    select_if(~ !any(is.na(.))) 
  
  # output results ------
  if (e_level == "METAB") { 
    OUT.FILE <- file.path(OUT_DIR, paste("lt", temp_feature, "temp_expr.tsv", sep = "_"))
  } else { 
    OUT.FILE <- file.path(OUT_DIR_PATHWAY, paste("lt_PW", st, temp_feature, "temp_expr.tsv", sep = "_"))
  }
  readr::write_tsv(metab_temp, file = OUT.FILE)
}



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()


