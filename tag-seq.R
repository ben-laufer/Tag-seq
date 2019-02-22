#3’ Tag RNA-seq pipeline

#References:
#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
#https://www.bioconductor.org/packages/release/bioc/vignettes/biobroom/inst/doc/biobroom_vignette.html
#https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html

# Load packages -----------------------------------------------------------

setwd("/Users/blaufer/Box Sync/PEBBLES")

rm(list=ls())
options(scipen=999)

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
#BiocManager::install(c("edgeR", "biobroom", "tidyverse", "stephenturner/annotables", "gplots", "RColorBrewer", "enrichR", "openxlsx"))
stopifnot(suppressMessages(sapply(c("edgeR", "biobroom", "tidyverse", "annotables", "gplots", "RColorBrewer", "enrichR", "openxlsx"), require, character.only = TRUE)))

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
  magrittr::set_rownames(ensemblIDs)
  
countMatrix <- countMatrix[-c(1:4),]

# Design Matrix -----------------------------------------------------------

designMatrix <- read.csv("sample_info.csv") %>%
  as_tibble()

# Tidy
designMatrix$Name <- gsub("-", "_", designMatrix$Name)
designMatrix$Name <- gsub(" ", "_", designMatrix$Name)
designMatrix$Name <- gsub("\n", "", designMatrix$Name)

# Treatment as numeric
designMatrix$Treatment <- as.character(designMatrix$Treatment)
designMatrix$Treatment[designMatrix$Treatment  == "Control"] <- "0"
designMatrix$Treatment <- as.numeric(designMatrix$Treatment)

# # Recode sex
# designMatrix$Sex <- as.character(designMatrix$Sex)
# designMatrix$Sex[designMatrix$Sex  == "F"] <- "0"
# designMatrix$Sex[designMatrix$Sex  == "M"] <- "1"
# designMatrix$Sex <- as.factor(designMatrix$Sex)

samples.idx <- pmatch(designMatrix$Name, colnames(countMatrix))
designMatrix <- designMatrix[order(samples.idx),]

# Preprocessing -----------------------------------------------------------

countMatrix <- countMatrix %>%
  dplyr::select(contains("P"))

designMatrix <- designMatrix %>%
  dplyr::filter(Tissue == "Placenta")

# Calculate normalization factors
countMatrix <- countMatrix %>%
  DGEList() %>%
  calcNormFactors()

# Filter genes with low expression

# # Automated method from ref 3
# keep.exprs <- filterByExpr(countMatrix, group = designMatrix$Treatment, lib.size = countMatrix$samples$lib.size)
# countMatrix <- countMatrix[keep.exprs,, keep.lib.sizes=FALSE]
# dim(countMatrix)

# Blythe method
cutoff <- 1
drop <- which(apply(cpm(countMatrix), 1, max) < cutoff)
countMatrix <- countMatrix[-drop,]
dim(countMatrix)

# Reorder design matrix 
samples.idx <- pmatch(designMatrix$Name, rownames(countMatrix$samples))
designMatrix <- designMatrix[order(samples.idx),]
stopifnot(rownames(countMatrix$samples) == designMatrix$Name)

# MDS of all interactions
group <- interaction(designMatrix$Treatment, designMatrix$Sex)
plotMDS(countMatrix, col = as.numeric(group))

# MDS of treatments simplified
plotMDS(countMatrix, col = (as.numeric(designMatrix$Treatment == "Control")+1))

# Voom transformation and calculation of variance weights -----------------

#mm <- model.matrix(~0 + designMatrix$Treatment + designMatrix$Sex) # Force zero intercept?
mm <- model.matrix(~designMatrix$Treatment + designMatrix$Sex)
logCPM <- voom(countMatrix, mm, plot = T)

# Fitting linear models in limma ------------------------------------------

fit <- lmFit(logCPM, mm)
head(coef(fit))

DEGs <- fit %>%
  contrasts.fit(coef = 2) %>%
  eBayes() %>%
  topTable(sort.by = "P", n = Inf) %>%
  rownames_to_column() %>% 
  as_tibble() %>%
  dplyr::rename(ensembl = rowname) %>% 
  dplyr::inner_join(grcm38, by = c("ensembl" = "ensgene")) %>% 
  dplyr::select(symbol, logFC, P.Value, adj.P.Val, ensembl, description) %>%
  dplyr::filter(P.Value < 0.05)
  
DEGs

# Plots -------------------------------------------------------------------

# Heatmap

heatMatrix <- logCPM$E[which(rownames(logCPM$E) %in% DEGs$ensembl),]

pdf("heatmap.pdf", height = 8.5, width = 11)

heatmap.2(heatMatrix,
          scale = "row",
          labCol = designMatrix$Treatment, 
          col = rev(brewer.pal(11, name = "RdBu")),
          trace = "none",
          dendrogram = "column")

dev.off()

# MD plot

dt <- fit %>% decideTests() %>% summary()

pdf("MDplot.pdf", height = 8.5, width = 11)

plotMD(fit,
       column = 2,
       status = dt[,2],
       main = colnames(fit)[2])

dev.off()

# Ontologies and Pathways -------------------------------------------------

# Check available databases
#dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018",
         "GO_Cellular_Component_2018",
         "GO_Molecular_Function_2018",
         "KEGG_2016",
         "Panther_2016",
         "Reactome_2016",
         "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")

GO <- DEGs %>%
  dplyr::select(symbol) %>%
  purrr::flatten() %>%
  enrichr(dbs)

write.xlsx(GO, file = "enrichr.xlsx", sep="")
