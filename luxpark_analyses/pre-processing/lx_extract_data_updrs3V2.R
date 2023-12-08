# Title: lx_extract_data_updrs3V2.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script filters all necessary data of PD-only subjects having UPDRS3 scores at V0, V1, V2.
# Usage: Rscript lx_extract_data_updrs3V2.R 
# Data: data from expression counts at timepoints V0, V1,V2, clinical targets data (UPDRS3) at V2.

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
analysis_name <- "02-pred-TS-UPDRS3-class"
IN_DIR <- "../data/00-cleansing/"
IN_DIR_PATHWAY <- "../data/01-dea-TS-PD/04-pathway_level"
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 
PHENO.FILE <- file.path(IN_DIR, "lx_clinical_targets_PD.tsv")
PHENO.OUT.FILE <- file.path(OUT_DIR, "lx_pheno.tsv") # PD at baseline, UPDRS3 (and relevant other vars) at multiple timepoints



# Main -------------------------------------------------------------------------
if ((!dir.exists(OUT_DIR)) | (!dir.exists(OUT_DIR_PATHWAY))) {
  dir.create(OUT_DIR, recursive = T)
  dir.create(OUT_DIR_PATHWAY, recursive = T)
}



# Data transformation ----------------------------------------------------------
pheno <- vroom(PHENO.FILE, show_col_types = FALSE) 

# filter PD patients + having data for V0, V1, V2 + having UPDRS3 info for all 3 timepoints
patient_list <-pheno %>%
  filter (DIAGNOSIS == "PD" & VISIT %in% c("V0", "V1", "V2")) %>%
  group_by(PATIENT_ID) %>% 
  dplyr::summarise(n = n(),
            na_updrs3 = sum(is.na(UPDRS__3))) %>%
  filter(n == 3 & na_updrs3 == 0) %>%
  pull(PATIENT_ID)

# filter PD patients from expression data
all_files <- c(list.files(path=IN_DIR, pattern = "log_transformed_data_fte.tsv", full.names=TRUE),
               list.files(path=IN_DIR_PATHWAY, pattern = "log_transformed_PW.*\\.tsv$", full.names=TRUE))

for (f in all_files) {
  filename <- str_extract(f, '[^/]*$')
  if (!str_detect(filename, pattern = "PW")) {
    OUT.FILE <- file.path(OUT_DIR, filename)
  } else {
    OUT.FILE <- file.path(OUT_DIR_PATHWAY, filename)
  }
  data <- vroom(f, show_col_types = FALSE) 
  data <- data %>% 
    dplyr::filter(PATIENT_ID %in% patient_list) 

  # output results ------
  readr::write_tsv(data, file = OUT.FILE) 
}

pheno <- pheno %>% 
  filter(PATIENT_ID %in% patient_list) %>% # patients with 3 timepoints
  filter (VISIT == "V2") # keep only phenotypical (i.e., target) features of V2
readr::write_tsv(pheno, file = PHENO.OUT.FILE)   



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()


