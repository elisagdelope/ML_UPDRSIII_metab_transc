# Title: ppmi_extract_data_updrs3v08.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script filters all necessary data of PD-only subjects having UPDRS3 score at V08 
# Usage: Rscript ppmi_extract_data_updrs3v08.R
# Data: data from temporal expression & phenotype from a previous analysis (02-pred-TS-PD)

# Packages ---------------------------------------------------------------------
library(readr)
library(dplyr)
library(vroom)
library(tidyr)
library(tibble)
library(stringr)



# I/O --------------------------------------------------------------------------
analysis_name <- "02-pred-TS-UPDRS3-class" # 02-pred-TS-UPDRS3 for regression
IN_DIR <- paste0("../data/02-pred-TS-PD/02-outfiles/temporal_data/") 
IN_DIR_PATHWAY <- paste0("../data/02-pred-TS-PD/04-pathway_level/temporal_data/") 
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles/temporal_data") 
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level/temporal_data") 
PHENO.FILE <- "../data/01-dea-PDHC/02-outfiles/ppmi_pheno_dsq.tsv"
PHENO.OUT.FILE <- paste0("../data/", analysis_name , "/02-outfiles/ppmi_pheno.tsv")



# Main -------------------------------------------------------------------------
if ((!dir.exists(OUT_DIR)) | (!dir.exists(OUT_DIR_PATHWAY))) {
  dir.create(OUT_DIR, recursive = T)
  dir.create(OUT_DIR_PATHWAY, recursive = T)
}



# Data transformation ----------------------------------------------------------
pheno <- vroom(PHENO.FILE, col_types = c("ccffdfildddddddddl")) 

# filter PD patients from pheno data with UPDRS3 at V08
PD_updrs3v08_list <- pheno %>% 
  rename_with(toupper) %>%
  filter (DIAGNOSIS == "PD" & VISIT == "V08" & !is.na(UPDRS3)) %>%
  pull(GSEID)

# filter PD patients from expression data
all_files <- c(list.files(path=IN_DIR, full.names=TRUE),
               list.files(path=IN_DIR_PATHWAY, full.names=TRUE))
var_id = "PATIENT_ID"

for (f in all_files) {
  
  data <- vroom(f, show_col_types = FALSE) 
  data <- data %>%
    filter (PATIENT_ID %in% PD_updrs3v08_list)

  # output results ------
  filename <- str_extract(f, '(?<=//).*')
  if (str_detect(filename, pattern = "GENE")) {
    OUT.FILE <- file.path(OUT_DIR, filename)
  } else {
    OUT.FILE <- file.path(OUT_DIR_PATHWAY, filename)
  }
  readr::write_tsv(data, file = OUT.FILE)   
}

pheno <- pheno %>%
  rename_with(toupper) %>%
  rename("GSEID" = "PATIENT_ID") %>% 
  filter(PATIENT_ID %in% PD_updrs3v08_list) %>% # patients with updrs3 at v08 time
  filter (VISIT == "V08") %>%
  filter(PATIENT_ID %in% data$PATIENT_ID)  # patients with temporal information

readr::write_tsv(pheno, file = PHENO.OUT.FILE)   



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()




