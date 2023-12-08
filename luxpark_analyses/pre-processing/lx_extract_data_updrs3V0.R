# Title: lx_extract_data_updrs3V0.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script filters all necessary data of PD-only subjects having UPDRS3 score at V0.
# Usage: Rscript lx_extract_data_updrs3V0.R 
# Data: data from expression counts at specific timepoint (V0), clinical targets data (UPDRS3) at V0.

# GC ---------------------------------------------------------------------------
rm(list = ls())
gc(T)



# Packages ---------------------------------------------------------------------
library(readr)
library(dplyr)
library(vroom)
library(tidyr)
library(tibble)
library(stringr)



# I/O --------------------------------------------------------------------------
analysis_name <- "02-pred-V0-UPDRS3-class"
IN_DIR <- paste0("../data/02-pred-V0-PD/02-outfiles/") 
IN_DIR_PATHWAY <- paste0("../data/02-pred-V0-PD/04-pathway_level/") 
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 

PHENO.FILE <- "../data/00-cleansing/lx_clinical_targets_PD.tsv"
PHENO.OUT.FILE <- file.path(OUT_DIR, "lx_pheno.tsv") # PD at baseline, UPDRS3 (and relevant other vars) at BL



# Main -------------------------------------------------------------------------
if ((!dir.exists(OUT_DIR)) | (!dir.exists(OUT_DIR_PATHWAY))) {
  dir.create(OUT_DIR, recursive = T)
  dir.create(OUT_DIR_PATHWAY, recursive = T)
}



# Data transformation ----------------------------------------------------------
pheno <- vroom(PHENO.FILE, show_col_types = FALSE) 

# filter PD patients from pheno data with UPDRS3 at V0
PD_updrs3V0_list <- pheno %>%
  dplyr::filter(VISIT == "V0" & !is.na(UPDRS__3)) %>%
  pull(SAMPLE_ID)

# filter PD patients from expression data
all_files <- c(list.files(path=IN_DIR, pattern = "log_transformed.*\\.tsv$", full.names=TRUE),
               list.files(path=IN_DIR_PATHWAY, pattern = "log_transformed_PW.*\\.tsv$", full.names=TRUE))

for (f in all_files) {
  filename <- str_extract(f, '(?<=//).*')
  if (!str_detect(filename, pattern = "PW")) {
    OUT.FILE <- file.path(OUT_DIR, filename)
  } else {
    OUT.FILE <- file.path(OUT_DIR_PATHWAY, filename)
  }
  data <- vroom(f, show_col_types = FALSE) 
  data <- data %>% 
    dplyr::filter(SAMPLE_ID %in% PD_updrs3V0_list) 

  # output results ------
  readr::write_tsv(data, file = OUT.FILE) 
}

pheno <- pheno %>% 
  filter(SAMPLE_ID %in% PD_updrs3V0_list)
readr::write_tsv(pheno, file = PHENO.OUT.FILE)   



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()


