## Load packages 
BiocManager::install("clusterProfiler",force = T)
BiocManager::install("pathview")
BiocManager::install("enrichplot",orce = T)
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(readr)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(WGCNA)
library(clusterProfiler)
library(enrichplot)
library(organism, character.only = TRUE)
library(tidyverse)
library(pheatmap)
library(ggforce)

# create directories for plots
dir.create(file.path(getwd(),'plots'), showWarnings = FALSE)
dir.create(file.path(getwd(),'data_output'), showWarnings = FALSE)

## Load data
protein_data <- read_csv("data/NULISA_TF_data.csv")
protein_data = protein_data[,2:ncol(protein_data)]
clinical_data <- read_csv("data/NULISA_TF_clinicaldata.csv")
clinical_data = clinical_data[,2:ncol(clinical_data)]
DE_ALS_CTR <- read_csv("data/NULISA_TF_ALSvsCTRL_DE.csv")
DE_basal_reflex <- read_csv("data/NULISA_TF_basalvsreflex_DE.csv")

colnames(DE_ALS_CTR) = colnames(DE_basal_reflex) = c("protein","log2FC","AvExp","T-stat","p-value","adj p-value","B")