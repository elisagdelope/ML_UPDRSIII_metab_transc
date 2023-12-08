# Title: lx_friedman_UPDRS3.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script performs comparisons among model performance, data types performance, and feature selection methods.
# Usage: Rscript lx_friedman_UPDRS3.R 
# Usage (can be run from command line): Rscript lx_friedman_UPDRS3.R > ../reports/PRED-V0-TS-PD_results/stdout_friedman.txt 2>err/stderr_friedmans.txt
# Data: data from summary of results xlsx file.


# GC ---------------------------------------------------------------------------
rm(list = ls())
gc(T)


# Packages ---------------------------------------------------------------------
library(readr)
library(plyr)
library(dplyr)
library(vroom)
library(tidyr)
library(tibble)
library(stringr)
library(PMCMRplus)
library(readxl)


#### CHANGE FOR UPDRS3
# I/O --------------------------------------------------------------------------
IN_DIR <- "../reports/PRED-V0-TS-UPDRS_3_results"
IN_FILE <- file.path(IN_DIR, "summary_AUC_scores_BLTS.tsv")
CONFOUNDERS_FILE <- file.path(IN_DIR, "CONFOUNDERS_summary_results.csv")
ORIG_FILE <- file.path(IN_DIR, "summary_results_UPDRS__3_binary.xlsx")
OUT_DIR <- IN_DIR
OUT_DIR_PLOTS <- "../reports/PRED-V0-TS-UPDRS_3_results/plots"
source("func_friedman_connover_holm_nxn.R")
source("func_friedman.R")


# Data transformations ---------------------------------------------------------

if (file.exists(IN_FILE)){
  auc <- vroom(IN_FILE, show_col_types = FALSE)
} else {
  # Data transformations ---------------------------------------------------------
  data <- list()
  data[["bl"]] <- read_excel(ORIG_FILE, sheet = 1)
  data[["lm_time"]] <- read_excel(ORIG_FILE, sheet = 2)
  data[["sd"]] <- read_excel(ORIG_FILE, sheet = 3)
  
  for (ft in names(data)) {
    data[[ft]] <- data[[ft]] %>%
      filter(.[[1]] != "de novo") %>%
      rename_at(1,~"data_type") %>%
      dplyr::select(-2) %>%
      mutate(across(-1, as.double)) %>% 
      mutate(linearSVM=as.double(linearSVM))
  }
  
  an_datatype <- list()
  an_model <- list()
  type_model_auc <- list()
  for (i in names(data)) {
    an_datatype[[i]] <- as.data.frame(apply(data[[i]] %>% 
                                              column_to_rownames(var = "data_type"), 1, median))
    colnames(an_datatype[[i]]) <- "AUC"
    an_datatype[[i]] <- an_datatype[[i]] %>% 
      rownames_to_column(var = "data_type")
    
    an_model[[i]] <- as.data.frame(apply(data[[i]] %>% 
                                           column_to_rownames(var = "data_type") %>% t, 1, median))
    colnames(an_model[[i]]) <- "AUC"
    an_model[[i]] <- an_model[[i]] %>% 
      rownames_to_column(var = "model")
    
    type_model_auc[[i]] <- data[[i]] %>% 
      pivot_longer(cols = -data_type, names_to = "model", values_to = "AUC")
    
  }
  
  
  
  # export AUC scores table ------------------------------------------------------
  #  merge all AUC scores into full AUC table
  auc <- bind_rows(type_model_auc, .id = "fttype")
  # print full AUC table
  readr::write_tsv(auc, file = file.path(OUT_DIR, "summary_AUC_scores_BLTS.tsv"))
}


# Pairwise comparison within each feature type: models x models -----------------------------------------
print("################## MODELS COMPARISON ANALYSIS ##################")

for (ft in unique(auc$fttype)){
  print(paste("#---", ft, "---#"))
  type_model_auc_i <- auc %>%
    dplyr::filter(fttype==ft)
  frd <- friedman(groups=type_model_auc_i$model, 
                  blocks=type_model_auc_i$data_type,
                  score=type_model_auc_i$AUC)
  if (frd$p.value < 0.05) {
    fr_nxn <- friedman_BergmannHommel_nxn(groups=type_model_auc_i$model,  
                                          blocks=type_model_auc_i$data_type,
                                          score=type_model_auc_i$AUC)
    print(paste("Friedman + Bergmann-Hommel post-hoc test: ", names(which.min(fr_nxn$ranks)), 
                "performs significantly (PAV < 0.05) better than", paste(names(which(fr_nxn$test[names(which.min(fr_nxn$ranks)),] < 0.05)), collapse=", ")))
    # print plot
    pdf(file = file.path(OUT_DIR_PLOTS, paste0("Friedman_nxn_models_",ft,".pdf")), width = 11, height = 4)
    print(fr_nxn$plot)
    dev.off()
  } else{
    print("Friedman hypothesis test was not rejected")
  }
}


# Pairwise comparison: data types x data types ---------------------------------
print("################## FEATURE TYPES COMPARISON ANALYSIS ##################")
for (ft in unique(auc$fttype)){
  print(paste("#---", ft, "---#"))
  type_model_auc_i <- auc %>%
  dplyr::filter(fttype==ft)
  frd <- friedman(groups=type_model_auc_i$data_type,
                  blocks=type_model_auc_i$model,
                  score=type_model_auc_i$AUC)
  if (frd$p.value < 0.05) {
    fr_nxn <- friedman_BergmannHommel_nxn(groups=type_model_auc_i$data_type,
                                          blocks=type_model_auc_i$model,
                                          score=type_model_auc_i$AUC)
    
    print(paste("Friedman + Bergmann-Hommel post-hoc test: ", names(which.min(fr_nxn$ranks)), 
                "performs significantly (PAV < 0.05) better than", paste(names(which(fr_nxn$test[names(which.min(fr_nxn$ranks)),] < 0.05)), collapse=", ")))
    # print plot
    pdf(file = file.path(OUT_DIR_PLOTS, paste0("Friedman_nxn_datatypes_",ft,".pdf")), width = 7, height = 4)
    print(fr_nxn$plot)
    dev.off()
    
  } else{
    print("Friedman hypothesis test was not rejected")
  }
}



# BL vs TS for genes -----------------------------------------------------------
print("################## BL VS TS IN METABOLITES COMPARISON ANALYSIS ##################")
auc_BLTS <- auc

type_model_auc_metabs <- auc_BLTS %>%
  dplyr::filter(data_type=="metabolites") %>%
  mutate(fttype = toupper(fttype))

print(paste("feature types compared:", paste(unique(type_model_auc_metabs$fttype), collapse=", ")))

frd <- friedman(groups=type_model_auc_metabs$fttype, 
                blocks=type_model_auc_metabs$model,
                score=type_model_auc_metabs$AUC)
if (frd$p.value < 0.05) {
  # LM (time) against the other features.
  #fr_1xn <- friedman_holm_nx1(groups=type_model_auc_metabs$fttype,  
  #                            blocks=type_model_auc_metabs$model,
  #                            score=type_model_auc_metabs$AUC, col_control=2)
  
  fr_nxn <- friedman_BergmannHommel_nxn(groups=type_model_auc_metabs$fttype,  
                                        blocks=type_model_auc_metabs$model,
                                        score=type_model_auc_metabs$AUC)
  print(paste("Friedman + Bergmann-Hommel post-hoc test: ", names(which.min(fr_nxn$ranks)), 
              "performs significantly (PAV < 0.05) better than", paste(names(which(fr_nxn$test[names(which.min(fr_nxn$ranks)),] < 0.05)), collapse=", ")))
  # print plot
  fr_nxn$plot +  
    scale_x_discrete(labels = function(labels) toupper(labels)) +
    scale_y_continuous(labels = function(labels) toupper(labels))
  pdf(file = file.path(OUT_DIR_PLOTS, "Friedman_nxn_BLTS_metabolitees.pdf"), width = 7, height = 4)
  print(fr_nxn$plot)
  dev.off()
} else{
  print("Friedman hypothesis test was not rejected")
  
}



# BL vs TS across data types (best model/data type) ----------------------------
print("################## BL VS TS across data types (best model/data type) COMPARISON ANALYSIS ##################")

type_model_auc_blts <- auc_BLTS %>%
  group_by(data_type, fttype) %>%
  slice_max(AUC) %>%
  mutate(fttype = toupper(fttype))

print(paste("feature types compared:", paste(unique(type_model_auc_blts$fttype), collapse=", ")))

frd <- friedman(groups=type_model_auc_blts$fttype, 
                blocks=type_model_auc_blts$data_type,
                score=type_model_auc_blts$AUC)
if (frd$p.value < 0.05) {
  fr_nxn <- friedman_BergmannHommel_nxn(groups=type_model_auc_blts$fttype,  
                                        blocks=type_model_auc_blts$data_type,
                                        score=type_model_auc_blts$AUC)
  print(paste("Friedman + Bergmann-Hommel post-hoc test: ", names(which.min(fr_nxn$ranks)), 
              "performs significantly (PAV < 0.05) better than", paste(names(which(fr_nxn$test[names(which.min(fr_nxn$ranks)),] < 0.05)), collapse=", ")))
  # print plot
  fr_nxn$plot
  pdf(file = file.path(OUT_DIR_PLOTS, "Friedman_nxn_BLTS_datatypes.pdf"), width = 7, height = 4)
  print(fr_nxn$plot)
  dev.off()
} else{
  print("Friedman hypothesis test was not rejected")
  
}



# within feature type (e.g. BL), comparison of DB-st RANKS ---------------------

frd <- friedman(groups=type_model_auc_blts$data_type, 
                blocks=type_model_auc_blts$fttype,
                score=type_model_auc_blts$AUC)
RM_TS <- list()
for (tf in c("BL", "LM_TIME", "SD")) {
  type_model_auc_blts_r <- type_model_auc_blts %>%
    dplyr::filter(fttype==tf) %>%
    dplyr::select(-model) %>%
    pivot_wider(names_from=data_type, values_from = AUC)
  rm <- colMeans(scmamp::rankMatrix(type_model_auc_blts_r[,-1]))
  rm <- data.frame(rm)
  colnames(rm) <- tf
  RM_TS[[tf]] <- rm
}
rank_df <- do.call(cbind, RM_TS)
rank_df$AVG_RANK <- rowMeans(rank_df)
rank_df <- rank_df[order(rank_df$AVG_RANK), ]


# plot ranks
#df <- reshape2::melt(fr_1xn$test[2,])
df <- reshape2::melt(rank_df %>%
                       rownames_to_column("data_type") %>%
                       pivot_longer(cols= c("BL", "LM_TIME", "SD", "AVG_RANK"))) %>%
  select(-variable)
colnames(df) <- c("data_type", "fttype", "value")

l <- rownames(rank_df)[order(rank_df$AVG_RANK, decreasing = TRUE)]
df$data_type <- factor(df$data_type, levels = l)
df <- df[order(df$data_type), ]

rank_plot <- ggplot2::ggplot(df, ggplot2::aes(x = factor(fttype, levels = c("AVG_RANK", "BL", "LM_TIME", "SD", "LM_LAG")), y = data_type, 
                                              fill = value)) + ggplot2::geom_tile(col = "white") + 
  ggplot2::scale_fill_continuous("value") + ggplot2::labs(x = "Algorithm", 
                                                          y = "Algorithm") + 
  ggplot2::geom_text(ggplot2::aes(label = round(value, 2)), size = 5, col = "white") + 
  labs(title="Ranked performance of database-pooling operator\ncombinations across baseline and temporal features") + 
  scale_fill_gradientn("Rank" , colours = c("#f9dc24", "#b73779", "#0d0887")) + # c("orange" , "skyblue4", "black") 
  labs(x = "Feature type", y = "Data type") + theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),    # Set x-axis label size
    axis.title.y = element_text(size = 14),     # Set y-axis label size
    plot.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()) 

# print plot
pdf(file = file.path(OUT_DIR_PLOTS, paste0("Friedman_ranks_aggs_BLTS.pdf")), width = 7, height = 4.5)
print(rank_plot)
dev.off()



# plot 1xn GOBP_mean against all rest---------
fr_1xn <- friedman_holm_nx1(groups=type_model_auc_blts$data_type,
                            blocks=type_model_auc_blts$fttype,
                            score=type_model_auc_blts$AUC,
                            col_control=3)

df <- reshape2::melt(fr_1xn$test[2,])
l <- colnames(fr_1xn$test)[order(fr_1xn$ranks, decreasing = TRUE)]
df$variable <- factor(df$variable, levels = l)
df <- df[order(df$variable), ]
df$Y <- df$variable[which(is.na(df$value))]
df <- df %>%
  filter(variable !=df$variable[which(is.na(df$value))])
plot <- ggplot2::ggplot(df, ggplot2::aes(x = variable, y = Y, 
                                         fill = value)) + ggplot2::geom_tile(col = "white") + 
  ggplot2::scale_fill_continuous("value") + ggplot2::labs(x = "Algorithm", 
                                                          y = "Algorithm") + 
  ggplot2::geom_text(ggplot2::aes(label = round(value, 2)), size = 5, col = "white") + 
  labs(title="") + 
  scale_fill_gradientn("Adj. PV" , colours = c("#f9dc24", "#b73779", "#0d0887")) + # c("orange" , "skyblue4", "black") 
  labs(x = "Data type", y = "") + theme(
    axis.text.x = element_text(size = 12, angle = 14, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),    # Set x-axis label size
    axis.title.y = element_text(size = 14),     # Set y-axis label size
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.title = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()) 


# print plot
pdf(file = file.path(OUT_DIR_PLOTS, paste0("Friedman_1xn_GOBP_mean.pdf")), width = 14, height = 2)
print(plot)
dev.off()





