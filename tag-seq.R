#3’ Tag RNA-seq pipeline

#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
#https://www.bioconductor.org/packages/release/bioc/vignettes/biobroom/inst/doc/biobroom_vignette.html

# Load packages -----------------------------------------------------------

rm(list=ls())
options(scipen=999)

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
#BiocManager::install(c("edgeR", "biobroom"))
stopifnot(suppressMessages(sapply(c("edgeR", "biobroom", "tidyverse"), require, character.only = TRUE)))

# Load data ---------------------------------------------------------------

setwd("/Users/blaufer/Box Sync/PEBBLES")
meta <- read.csv("sample_info.csv")
files <- list.files(path = glue::glue(getwd(), "/GeneCounts"), pattern = "*.ReadsPerGene.out.tab")


# Format count matrix -----------------------------------------------------

name <- gsub( "(?:[^_]+_){4}([^_ ]+)*$","", list.files(glue::glue(getwd(), "/GeneCounts"), pattern = "*.ReadsPerGene.out.tab"))

