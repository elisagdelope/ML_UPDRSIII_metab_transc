# Title: ppmi_extract_data_updrs3BL.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script filters all necessary data of PD-only subjects having UPDRS3 score at BL 
# Usage: Rscript ppmi_extract_data_updrs3BL.R
# Data: data from temporal expression & phenotype from a previous analysis (02-pred-BL-PD)

# Packages ---------------------------------------------------------------------
library(readr)
library(dplyr)
library(vroom)
library(tidyr)
library(tibble)
library(stringr)



# I/O --------------------------------------------------------------------------
analysis_name <- "02-pred-BL-UPDRS3-class" # 02-pred-BL-UPDRS3 for regression
IN_DIR <- paste0("../data/02-pred-BL-PD/02-outfiles/") 
IN_DIR_PATHWAY <- paste0("../data/02-pred-BL-PD/04-pathway_level/") 
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 
PHENO.FILE <- "../data/01-dea-PDHC/02-outfiles/ppmi_pheno_dsq.tsv"
PHENO.OUT.FILE <- file.path(OUT_DIR, "ppmi_pheno.tsv")



# Main -------------------------------------------------------------------------
if ((!dir.exists(OUT_DIR)) | (!dir.exists(OUT_DIR_PATHWAY))) {
  dir.create(OUT_DIR, recursive = T)
  dir.create(OUT_DIR_PATHWAY, recursive = T)
}



# Data transformation ----------------------------------------------------------
pheno <- vroom(PHENO.FILE, col_types = c("ccffdfildddddddddl")) 

# filter PD patients from pheno data with UPDRS3 at BL
PD_updrs3vBL_list <- pheno %>% 
  rename_with(toupper) %>%
  filter (DIAGNOSIS == "PD" & VISIT == "BL" & !is.na(UPDRS3)) %>%
  pull(GSEID)

# filter PD patients from expression data
all_files <- c(list.files(path=IN_DIR, pattern = "flt_.*\\.tsv$", full.names=TRUE),
               list.files(path=IN_DIR_PATHWAY, pattern = "_expression.tsv$", full.names=TRUE))

for (f in all_files) {
  
  filename <- str_extract(f, '(?<=//).*')
  if (str_detect(filename, pattern = "star")) {
    OUT.FILE <- file.path(OUT_DIR, filename)
    var_id = "GENEID"
  } else {
    OUT.FILE <- file.path(OUT_DIR_PATHWAY, filename)
    var_id = str_extract(f, '(?<=//).+?(?=_)') 
    var_id = paste0(var_id, "_NAME") # retrieve "db_NAME" info (e.g. GOCC_NAME)
  }
  data <- vroom(f, show_col_types = FALSE) 
  data <- data %>% 
    rename_with(toupper) %>%
    column_to_rownames(var_id) %>%
    dplyr::select(starts_with(PD_updrs3vBL_list)) %>%
    rownames_to_column(var_id)
  
  # output results ------
  readr::write_tsv(data, file = OUT.FILE) 
}

pheno <- pheno %>% 
  filter(GSEID %in% PD_updrs3vBL_list & visit == "BL") %>%
  rename_with(toupper) %>%
  rename("GSEID" = "PATIENT_ID")
  
readr::write_tsv(pheno, file = PHENO.OUT.FILE)   



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()



