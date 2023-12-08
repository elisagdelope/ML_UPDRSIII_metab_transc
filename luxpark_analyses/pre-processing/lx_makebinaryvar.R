# Title: lx_makebinaryvar.R 
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script transforms UPDRS__3 continuous variable from continuous to binary (0,1) correspondingly to (below median, above median) and adds the new binary variable to the pheno file.
# Info: Pheno file is overwritten -> one time run only!
# Usage : Rscript lx_makebinaryvar.R 
# Data: clinical data with continuous variable to transform into binary.


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
IN_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
IN_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 
OUT_DIR <- IN_DIR
OUT_DIR_PATHWAY <- IN_DIR_PATHWAY
PHENO.FILE <- file.path(IN_DIR, "lx_pheno.tsv") # clinical vars, including target UPDRS__3
target <- "UPDRS__3"



# Main -------------------------------------------------------------------------
if ((!dir.exists(OUT_DIR)) | (!dir.exists(OUT_DIR_PATHWAY))) {
  dir.create(OUT_DIR, recursive = T)
  dir.create(OUT_DIR_PATHWAY, recursive = T)
}



# Data transformation ----------------------------------------------------------
pheno <- vroom(PHENO.FILE, show_col_types = FALSE) 
threshold <- median(pheno[[target]])
binary_target <- paste0(target, "_binary")
#pheno <- pheno %>%
#  mutate(!!sym(binary_target) := ifelse(!!sym(target) < threshold, 0, 1)) 
if (sum(pheno[[target]] > threshold) > sum(pheno[[target]] < threshold)) {
  pheno[[binary_target]] <- ifelse(pheno[[target]] <= threshold, 0, 1) # samples with threshold go to the negative class if upper segment is more frequent
} else {
  pheno[[binary_target]] <- ifelse(pheno[[target]] < threshold, 0, 1) # samples with threshold go to the positive class if lower segment is more frequent
}


# export outputs  --------------------------------------------------------------
readr::write_tsv(pheno, file = PHENO.FILE) # over-write pheno file with new variable



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()


