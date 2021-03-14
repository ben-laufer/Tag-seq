#!/usr/bin/env Rscript

# WGCNA
# Modified version of WGCNA tutorials with custom functions
# Ref: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html

# Install and Update ------------------------------------------------------

rm(list=ls())
options(scipen=999)
options(stringsAsFactors = FALSE)

cat("\n[WGCNA] Installing and updating pacakges \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

CRAN <- c("BiocManager", "remotes")
new.packages <- CRAN[!(CRAN %in% installed.packages()[,"Package"])]
if(length(new.packages)>0){
  install.packages(new.packages, repos ="https://cloud.r-project.org")
}
stopifnot(suppressMessages(sapply(CRAN, require, character.only= TRUE)))

packages <- c("matrixStats","Hmisc","splines","foreach","doParallel","fastcluster","dynamicTreeCut","survival","GO.db","preprocessCore","impute","WGCNA","tidyverse")

new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0){
  library("BiocManager")
  BiocManager::install(new.packages, ask = FALSE)
}
stopifnot(suppressMessages(sapply(packages, require, character.only= TRUE)))
suppressWarnings(BiocManager::valid(fix = TRUE, update = TRUE, ask = FALSE))

# Global variables --------------------------------------------------------

cores <- 30
# WGCNA specific options
options(stringsAsFactors = FALSE)
enableWGCNAThreads(cores)

## WGCNA 1: Data input and cleaning ----
setwd("/Users/blaufer/Box Sync/PEBBLES/tag-seq")

logCPM <- openxlsx::read.xlsx("Placenta_voomLogCPMforWGCNA.xlsx")
#logCPM <- openxlsx::read.xlsx("Brain_voomLogCPMforWGCNA.xlsx") 

# Tidy
WGCNA_data0 <- logCPM %>%
  remove_rownames() %>% 
  column_to_rownames("Gene") %>% 
  as.matrix() %>% 
  t()

# Check for missing values and outliers
WGCNA_gsg <- goodSamplesGenes(WGCNA_data0, verbose = 3)
WGCNA_gsg$allOK

# Remove any offending regions and/or samples
if (!WGCNA_gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!WGCNA_gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(WGCNA_data0)[!WGCNA_gsg$goodGenes], collapse = ", ")));
  if (sum(!WGCNA_gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(WGCNA_data0)[!WGCNA_gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  WGCNA_data0 = WGCNA_data0[WGCNA_gsg$goodSamples, WGCNA_gsg$goodGenes]
}

# Cluster samples to check for outliers
sampleTree = hclust(dist(WGCNA_data0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

## Manually remove outliers based on visualization (Change for each analysis)
# Plot a line to show what the cut would exclude
abline(h = 120, col = "red"); # 120 for placenta
# Determine cluster under the line (use same value as line to cut)
clust = cutreeStatic(sampleTree, cutHeight = 120, minSize = 10) #120 for placenta
table(clust)

# clust 1 contains the samples we want to keep, so remove outliers here (clust problem! minsize?)
keepSamples = (clust==1)
WGCNA_data = WGCNA_data0[keepSamples, ]
nGenes = ncol(WGCNA_data)
nSamples = nrow(WGCNA_data)

# Reload covariates (How do I do this automagically?)
#indx <- sapply(meta, is.factor)
#test <- lapply(meta[indx], function(x) as.numeric(levels(x))[x])
meta <- read.csv("sample_info.csv", header = TRUE)

# Tidy
meta$Name <- gsub("-", "_", meta$Name)
meta$Name <- gsub(" ", "_", meta$Name)
meta$Name <- gsub("\n", "", meta$Name)

# Treatment as a categorical (factor) since continous (numeric) assumes a linear response and PCB response is non-monotonic
meta$Treatment <- as.character(meta$Treatment)
meta$Treatment[meta$Treatment  == "Control"] <- "0"
meta$Treatment <- as.factor(meta$Treatment)

# Add litter
meta$Litter <- str_split_fixed(meta$Name, "_", n = 3)[,1] %>%
  as.factor()

# Recode sex
meta$Sex <- as.character(meta$Sex)
meta$Sex[meta$Sex  == "F"] <- "0"
meta$Sex[meta$Sex  == "M"] <- "1"
meta$Sex <- as.factor(meta$Sex)

# Filter for relevant samples
meta <- meta %>%
  filter(Tissue == "Brain") %>% # Placenta
  select(-Tissue)

meta <- meta[which(meta$Name %in% rownames(WGCNA_data)),]
samples.idx <- pmatch(meta$Name, rownames(WGCNA_data))
meta <- meta[order(samples.idx),]
str(meta)

# Error check for subsetting and order
stopifnot(rownames(WGCNA_data) == meta$Name)

WGCNA_samples <- rownames(WGCNA_data)
traitRows <- match(WGCNA_samples, meta$Name)
datTraits <- meta[traitRows, -1]

# Make factors numeric
datTraits$Treatment <- as.numeric(datTraits$Treatment)
datTraits$Sex <- as.numeric(datTraits$Sex)
datTraits$Litter <- as.numeric(datTraits$Litter)
datTraits <- datTraits[ , c(1,2,3)]

rownames(datTraits) <- meta[traitRows, 1]
collectGarbage()

# Re-cluster samples
sampleTree2 = hclust(dist(WGCNA_data), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = TRUE);

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")

dev.off()

## WGCNA 2: Automatic, one-step network construction and module detection ----

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(WGCNA_data, powerVector = powers, corFnc = "bicor", networkType = "unsigned", verbose = 5)

pdf("soft_thresholding_power.pdf", height = 5, width =9)
# Plot the results:
#sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

cat("\n[DM.R] Saving Rdata \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
WGCNA_1_env <- ls(all=TRUE)
save(list = WGCNA_1_env , file = "WGCNA_1.RData") 
#load("WGCNA_1.RData")


# Epigenerate -------------------------------------------------------------
# Run this part on epigenerate and request it for a day

cd /share/lasallelab/Ben/PEBBLES/tag-seq/placentaWGCNA or cd /share/lasallelab/Ben/PEBBLES/tag-seq/BrainWGCNA
module load R/3.6.0
aklog
R

load("WGCNA_1.RData")
library(WGCNA) # WGCNA needs to be the last loaded package https://www.biostars.org/p/305714/

#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

# Construct network and detect modules 
# Use identified soft threshold (power) above for power below or if no obvious plateau use Horvath's pre-caluclated (topic 6 from link below)
# https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html
net = blockwiseModules(WGCNA_data,
                       power = 5, # Choose from above plot, the lowest power where the curve flattens out upon reaching a high value (~0.90) # 5 for placenta, 7 for Brain
                       networkType = "unsigned",
                       TOMType = "unsigned", 
                       corType = "bicor", # More powerful than "pearson", https://www.ncbi.nlm.nih.gov/pubmed/23217028
                       maxPOutliers = 0.10, # Forces bicor to never regard more than the specified proportion of samples as outliers (https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html)
                       numericLabels = TRUE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "regionTOM", 
                       blocks = NULL,
                       maxBlockSize = length(WGCNA_data), # split calculations into blocks or don't using = length(WGCNA_data), limit is = sqrt(2^31)
                       nThreads = 25,
                       verbose = 3)

# open a graphics window
#sizeGrWindow(12, 9)
pdf("module_dendogram.pdf", height = 5, width =9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "WGCNA-networkConstruction-auto.RData")

# Back to the desktop -----------------------------------------------------

## WGCNA 3: Relating modules to phenotypes and identifying important regions ----
setwd("/Users/blaufer/Box Sync/PEBBLES/tag-seq/Placenta WGCNA")
#setwd("/Users/blaufer/Box Sync/PEBBLES/tag-seq/Brain WGCNA")

library(tidyverse)
library(WGCNA)

load("WGCNA_1.RData")
load("WGCNA-networkConstruction-auto.RData")
options(stringsAsFactors = FALSE)

# Quantify module-trait associations 
# Define numbers of genes and samples
nGenes = ncol(WGCNA_data)
nSamples = nrow(WGCNA_data)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(WGCNA_data, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


#' tidyWGCNA
#' @description remove NA values for factors with single levels WGCNA modules
#' @param moduleTrait Matrix of module trait p-values or correlation values from WGCNA tutorial 3
#' @return tidied module trait matrix
#' @export tidyWGCNA
tidyWGCNA <- function(moduleTrait = moduleTrait){
  moduleTrait <- moduleTrait  %>%
    as.data.frame()
  moduleTrait  <- moduleTrait[, unlist(lapply(moduleTrait, function(x) !all(is.na(x))))]
  moduleTrait <- moduleTrait %>%
    as.matrix %>%
    return()
}

moduleTraitPvalue <- tidyWGCNA(moduleTraitPvalue)
moduleTraitCor <- tidyWGCNA(moduleTraitCor)

#' fdrWGCNA
#' @description fdr correct correlations for significant WGCNA modules
#' @param moduleTraitPvalue Matrix of module trait p-values from WGCNA tutorial 3
#' @return fdr corrected module trait p-value matrix
#' @export fdrWGCNA
fdrWGCNA <- function(moduleTraitPvalue = moduleTraitPvalue){
  colNames <- colnames(moduleTraitPvalue)
  rowNames <- row.names(moduleTraitPvalue)
  moduleTraitPvalue <- moduleTraitPvalue %>% 
    as.matrix %>% 
    as.vector %>% 
    p.adjust(method='fdr') %>% 
    matrix(ncol=ncol(moduleTraitPvalue)) %>%
    as.data.frame()
  colnames(moduleTraitPvalue) <- colNames
  rownames(moduleTraitPvalue) <- rowNames
  moduleTraitPvalue <- moduleTraitPvalue %>%
    as.matrix() %>%
    return()
}

moduleTraitPvalue <- fdrWGCNA(moduleTraitPvalue)

# Remove zero variance from datTraits
datTraits <- datTraits %>% dplyr::select(one_of(colnames(moduleTraitCor)))

# Visualize relationships 
#sizeGrWindow(10,6)
pdf("module_trait_correlations.pdf", height = 6, width = 10)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


# DMR signficance and module membership
# Define variable diagnosis containing the diagnosis column of datTrait
diagnosis = as.data.frame(datTraits$Diagnosis);
names(diagnosis) = "diagnosis"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(WGCNA_data, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(WGCNA_data, diagnosis, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(diagnosis), sep="");
names(GSPvalue) = paste("p.GS.", names(diagnosis), sep="");

# Intramodular analysis
# Choose module from correlation heatmap
module = "blue"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for diagnosis",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# Create the starting data frame for saving
geneInfo0 = data.frame(moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for diagnosis
modOrder = order(-abs(cor(MEs, diagnosis, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.diagnosis))
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "WGCNA_DMR_Info.csv")

library(tidyverse)
library(GenomicRanges)

#' tidyModules
#' @description Extract and tidy regions of methylation from significant WGCNA modules
#' @param sigModules Character vector with names of signficant modules 
#' @param geneInfo geneInfo dataframe from WGCNA analysis (end of tutorial 3)
#' @return Genomic ranges object with coordinates for regions of methylation and module names as meta information
#' @export tidyModules
tidyModules <- function(sigModules = sigModules,
                        geneInfo = geneInfo){
  geneInfo %>%
    rownames_to_column() %>%
    as.tibble() %>%
    filter(moduleColor == sigModules) %>%
    dplyr::select(rowname, moduleColor) %>%
    separate(rowname, into = c('seqnames', 'start', 'end')) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE) %>%
    return()
}

sigModuleRanges <- tidyModules(sigModules = c("black", "yellow"), geneInfo = geneInfo)

sigModuleRanges <- split(sigModuleRanges, sigModuleRanges$moduleColor)

seqlevelsStyle(sigModuleRanges) <- "UCSC" 

black <- sigModuleRanges[[1]]
yellow <- sigModuleRanges[[2]]

backgroundRanges <- makeGRangesFromDataFrame(smooth)
seqlevelsStyle(backgroundRanges) <- "UCSC" 

library(liftOver)

path <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch <- import.chain(path)
ch

black_liftOver <- granges(unlist(liftOver(black, ch)))
yellow_liftOver <- granges(unlist(liftOver(yellow, ch)))
background_liftOver <- unlist(liftOver(backgroundRanges, ch))

library(rGREAT)
library(openxlsx)

black_GREAT <- submitGreatJob(black_liftOver,
                              bg = background_liftOver,
                              species = "hg19",
                              request_interval = 0)
black_tb <- getEnrichmentTables(black_GREAT, category = c("GO", "Pathway Data"))
write.xlsx(black_tb , file = "black_GREAT_results.xlsx", sep="")


yellow_GREAT <- submitGreatJob(yellow_liftOver,
                               bg = background_liftOver,
                               species = "hg19",
                               request_interval = 0)
yellow_tb <- getEnrichmentTables(yellow_GREAT, category = c("GO", "Pathway Data"))
write.xlsx(yellow_tb, file = "yellow_GREAT_results.xlsx", sep="")


WGNCA_env <- ls(all=TRUE)
save(list = WGNCA_env, file = "WGCNA.RData") 
#load("WGCNA.RData")

options(stringsAsFactors = TRUE)
