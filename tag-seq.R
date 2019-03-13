#3’ Tag RNA-seq pipeline

#References:
#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
#https://www.bioconductor.org/packages/release/bioc/vignettes/biobroom/inst/doc/biobroom_vignette.html
#https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html

# Load packages -----------------------------------------------------------

setwd("/Users/blaufer/Desktop/PEBBLES/tag-seq")

rm(list=ls())
options(scipen=999)

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
#BiocManager::install(c("edgeR", "biobroom", "tidyverse", "stephenturner/annotables", "gplots", "RColorBrewer", "enrichR", "openxlsx", "rstudio/gt", "plyr"))
stopifnot(suppressMessages(sapply(c("edgeR", "biobroom", "tidyverse", "annotables", "gplots", "RColorBrewer", "enrichR", "openxlsx", "gt", "plyr"), require, character.only = TRUE)))

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

# Could alternatively use edgeR::readDGE() but that calls to the slower read.delim()
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

# Add litter
designMatrix$Litter <- str_split_fixed(designMatrix$Name, "_", n =3)[,1] %>%
  as.factor()

# # Recode sex
# designMatrix$Sex <- as.character(designMatrix$Sex)
# designMatrix$Sex[designMatrix$Sex  == "F"] <- "0"
# designMatrix$Sex[designMatrix$Sex  == "M"] <- "1"
# designMatrix$Sex <- as.factor(designMatrix$Sex)

samples.idx <- pmatch(designMatrix$Name, colnames(countMatrix))
designMatrix <- designMatrix[order(samples.idx),]

# Preprocessing -----------------------------------------------------------

# Select sample subset
countMatrix <- countMatrix %>%
  dplyr::select(contains("P")) #B

designMatrix <- designMatrix %>%
  dplyr::filter(Tissue == "Placenta") #Brain

# Calculate normalization factors
countMatrix <- countMatrix %>%
  DGEList() %>%
  calcNormFactors()

# Raw density of log-CPM values

L <- mean(countMatrix$samples$lib.size) * 1e-6
M <- median(countMatrix$samples$lib.size) * 1e-6

logCPM <- cpm(countMatrix, log = TRUE)
logCPM.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(countMatrix)
col <- brewer.pal(nsamples, "Paired")

pdf("density_plot.pdf", height = 8.5, width = 11)
par(mfrow = c(1,2))
plot(density(logCPM[,1]), col = col[1], lwd = 2, las = 2, main = "", xlab = "")
title(main = "A. Raw data", xlab = "Log-cpm")
abline(v = logCPM.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(logCPM[,i])
  lines(den$x, den$y, col = col[i], lwd = 2)
}
legend("topright", designMatrix$Name, text.col = col, bty = "n", cex = 0.5)


# Filter genes with low expression

# Automated method from ref 3
keep.exprs <- filterByExpr(countMatrix, group = designMatrix$Treatment, lib.size = countMatrix$samples$lib.size)
countMatrix <- countMatrix[keep.exprs,, keep.lib.sizes=FALSE]
dim(countMatrix)

# # Blythe method for filtering
# cutoff <- 1
# drop <- which(apply(cpm(countMatrix), 1, max) < cutoff)
# countMatrix <- countMatrix[-drop,]
# dim(countMatrix)
# countMatrix$samples$lib.size <- colSums(countMatrix$counts) # Reset library size after filtering

# Reorder design matrix 
samples.idx <- pmatch(designMatrix$Name, rownames(countMatrix$samples))
designMatrix <- designMatrix[order(samples.idx),]
stopifnot(rownames(countMatrix$samples) == designMatrix$Name)

# Filtered density plot of log-CPM values 
logCPM <- cpm(countMatrix, log = TRUE)
plot(density(logCPM[,1]), col = col[1], lwd = 2, las =2 , main = "", xlab = "")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=logCPM.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(logCPM[,i])
  lines(den$x, den$y, col = col[i], lwd = 2)
}
legend("topright", designMatrix$Name, text.col = col, bty = "n", cex = 0.5)
dev.off()

# MDS of all interactions
group <- interaction(designMatrix$Treatment, designMatrix$Sex)
plotMDS(countMatrix, col = as.numeric(group))

# MDS of treatments simplified
plotMDS(countMatrix, col = (as.numeric(designMatrix$Treatment == "Control")+1))

# Voom transformation and calculation of variance weights -----------------

#mm <- model.matrix(~0 + designMatrix$Treatment + designMatrix$Sex) # Force zero intercept?
mm <- model.matrix(~designMatrix$Treatment + designMatrix$Sex +designMatrix$Litter)
voomLogCPM <- voom(countMatrix, mm, plot = T)

# Boxplots of logCPM values before and after normalization
par(mfrow=c(1,2))
boxplot(logCPM, las=2, col=col, main="")
title(main="A. Unnormalised data",ylab="Log-cpm")
boxplot(voomLogCPM$E, las=2, col=col, main="")
title(main="B. Normalised data",ylab="Log-cpm")

# Fitting linear models in limma ------------------------------------------

fit <- lmFit(voomLogCPM, mm)
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
  
# HTML report -------------------------------------------------------------

DEGs %>%
  dplyr::rename(Gene = symbol,
                "p-value" = P.Value,
                "q-value" = adj.P.Val,
                Description = description) %>% 
  dplyr::select(-ensembl) %>%
  dplyr::mutate(Description = purrr::map_chr(strsplit(DEGs$description, split='[', fixed=TRUE),function(x) (x[1]))) %>% 
  gt() %>%
  tab_header(
    title = glue::glue("{nrow(DEGs)} Differentially Expressed Genes"),
    subtitle = glue::glue("{round(sum(DEGs$logFC > 0) / nrow(DEGs), digits = 2)*100}% up-regulated, \\
                          {round(sum(DEGs$logFC < 0) / nrow(DEGs), digits = 2)*100}% down-regulated")
    ) %>% 
  fmt_number(
    columns = vars("logFC"),
    decimals = 2
  ) %>% 
  fmt_scientific(
    columns = vars("p-value", "q-value"),
    decimals = 2
  ) %>%
  as_raw_html(inline_css = TRUE) %>%
  write("DEGs.html") 

# Plots -------------------------------------------------------------------

# Heatmap

heatMatrix <- voomLogCPM$E[which(rownames(voomLogCPM$E) %in% DEGs$ensembl),]

# https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_color_hue <- function(n = n){
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

gg_color <- c(gg_color_hue(length(levels(as.factor(designMatrix$Treatment)))))

ColSideColors <- plyr::mapvalues(as.factor(designMatrix$Treatment),
                                 from = levels(as.factor(designMatrix$Treatment)),
                                 to = unique(gg_color)) 

pdf("heatmap.pdf", height = 8.5, width = 11)

heatmap.2(heatMatrix,
          scale = "row",
          labCol = designMatrix$Treatment,
          labRow = NA,
          col = rev(brewer.pal(11, name = "RdBu")),
          trace = "none",
          main = glue::glue("{nrow(DEGs)} Differentially Expressed Genes"),
          Rowv= as.dendrogram(hclust(dist(heatMatrix))),
          Colv = T,
          ColSideColors = as.character(ColSideColors))

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
