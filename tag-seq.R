#3’ Tag RNA-seq pipeline

#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
#https://www.bioconductor.org/packages/release/bioc/vignettes/biobroom/inst/doc/biobroom_vignette.html

# Load packages -----------------------------------------------------------

setwd("/Users/blaufer/Box Sync/PEBBLES")

rm(list=ls())
options(scipen=999)

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
#BiocManager::install(c("edgeR", "biobroom"))
stopifnot(suppressMessages(sapply(c("edgeR", "biobroom", "tidyverse"), require, character.only = TRUE)))

# Count Matrix ------------------------------------------------------------

#name <- gsub( "(?:[^_]+_){4}([^_ ]+)*$","", files)

# STAR quantMode geneCounts output:
#column 1: gene ID
#column 2: counts for unstranded RNA-seq
#column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
#column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)

sampleNames <- list.files(path = glue::glue(getwd(), "/GeneCounts"), pattern = "*.ReadsPerGene.out.tab") %>%
  stringr::str_split_fixed("_", n = 4) %>%
  as_tibble() %>%
  tidyr::unite(Name, c(V1:V3), sep = "_") %>%
  dplyr::select(Name) %>% 
  purrr::flatten_chr()

ensemblIDs <- list.files(path = glue::glue(getwd(), "/GeneCounts"), pattern = "*.ReadsPerGene.out.tab", full.names = T)[1] %>% 
  data.table::fread(select = 1) %>%
  purrr::flatten_chr()

countMatrix <- list.files(path = glue::glue(getwd(), "/GeneCounts"), pattern = "*.ReadsPerGene.out.tab", full.names = T) %>%
  purrr::map_dfc(data.table::fread, select = 3, data.table = F) %>%
  magrittr::set_colnames(sampleNames) %>% 
  magrittr::set_rownames(ensemblIDs) %>%
  as.matrix()
  
countMatrix <- countMatrix[-c(1:4),]

# Design Matrix -----------------------------------------------------------

designMatrix <- read.csv("sample_info.csv") %>%
  as_tibble()

designMatrix$Name <- gsub("-", "_", designMatrix$Name)
designMatrix$Name <- gsub(" ", "_", designMatrix$Name)
designMatrix$Name <- gsub("\n", "", designMatrix$Name)

# Preprocessing -----------------------------------------------------------

# Calculate normalization factors
countMatrix <- countMatrix %>%
  DGEList() %>%
  calcNormFactors()

# Filter genes with low expression
cutoff <- 1
drop <- which(apply(cpm(countMatrix), 1, max) < cutoff)
countMatrix <- countMatrix[-drop,] 
dim(countMatrix) 

# Group information
stopifnot(rownames(countMatrix$samples) == designMatrix$Name)

plotMDS(countMatrix)


