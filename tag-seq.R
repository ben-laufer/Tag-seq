#3’ Tag RNA-seq pipeline

#References:
#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
#https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html

# Load packages -----------------------------------------------------------

setwd("/Users/blaufer/Desktop/PEBBLES/tag-seq")

rm(list=ls())
options(scipen=999)

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
#BiocManager::install(c("edgeR", "tidyverse", "stephenturner/annotables", "gplots", "RColorBrewer", "enrichR", "openxlsx", "rstudio/gt", "plyr", "glue", "Glimma", "sva", "cowplot"))
stopifnot(suppressMessages(sapply(c("edgeR", "tidyverse", "annotables", "gplots", "RColorBrewer", "enrichR", "openxlsx", "gt", "plyr", "glue", "Glimma", "sva", "cowplot"),
                                  require, character.only = TRUE)))

for(tissue in 1:2){
  
  # Count Matrix ------------------------------------------------------------
  
  #name <- gsub( "(?:[^_]+_){4}([^_ ]+)*$","", files)
  
  # STAR quantMode geneCounts output:
  #column 1: gene ID
  #column 2: counts for unstranded RNA-seq
  #column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
  #column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)
  
  sampleNames <- list.files(path = glue::glue(getwd(), "/GeneCounts"), pattern = "*.ReadsPerGene.out.tab") %>%
    stringr::str_split_fixed("_", n = 4) %>%
    tibble::as_tibble() %>%
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
    tibble::as_tibble()
  
  # Tidy
  designMatrix$Name <- gsub("-", "_", designMatrix$Name)
  designMatrix$Name <- gsub(" ", "_", designMatrix$Name)
  designMatrix$Name <- gsub("\n", "", designMatrix$Name)
  
  # Treatment as a categorical (factor) since continous (numeric) assumes a linear response and PCB response is non-monotonic
  designMatrix$Treatment <- as.character(designMatrix$Treatment)
  designMatrix$Treatment[designMatrix$Treatment  == "Control"] <- "0"
  designMatrix$Treatment <- as.factor(designMatrix$Treatment)
  
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
  
  # Assign tissue variables
  tissueLetter <- case_when(tissue == 1 ~ "P",
                           tissue == 2 ~ "B")
  
  tissueName <- case_when(tissue == 1 ~ "Placenta",
                         tissue == 2 ~ "Brain")
  
  print(glue::glue("Preprocessing {tissueName} samples"))
  
  # Select sample subset
  countMatrix <- countMatrix %>%
    dplyr::select(contains(tissueLetter))
  
  designMatrix <- designMatrix %>%
    dplyr::filter(Tissue == tissueName)
  
  # Create DGE list and calculate normalization factors
  countMatrix <- countMatrix %>%
    DGEList() %>%
    calcNormFactors()
  
  # Reorder design matrix 
  samples.idx <- pmatch(designMatrix$Name, rownames(countMatrix$samples))
  designMatrix <- designMatrix[order(samples.idx),]
  stopifnot(rownames(countMatrix$samples) == designMatrix$Name)
  
  # Add sample info from design matrix to DGE list
  countMatrix$samples$group <- designMatrix$Treatment
  countMatrix$samples$Sex <- designMatrix$Sex
  countMatrix$samples$Litter <- designMatrix$Litter
  countMatrix$samples$Tissue <- designMatrix$Tissue
  
  # Raw density of log-CPM values
  
  L <- mean(countMatrix$samples$lib.size) * 1e-6
  M <- median(countMatrix$samples$lib.size) * 1e-6
  
  logCPM <- cpm(countMatrix, log = TRUE)
  logCPM.cutoff <- log2(10/M + 2/L)
  nsamples <- ncol(countMatrix)
  col <- brewer.pal(nsamples, "Paired")
  
  pdf(glue::glue("{tissueName}_density_plot.pdf"), height = 8.5, width = 11)
  par(mfrow = c(1,2))
  
  plot(density(logCPM[,1]), col = col[1], lwd = 2, las = 2, main = "", xlab = "")
  title(main = "A. Raw data", xlab = "Log-cpm")
  abline(v = logCPM.cutoff, lty = 3)
  for (i in 2:nsamples){
    den <- density(logCPM[,i])
    lines(den$x, den$y, col = col[i], lwd = 2)
  }
  legend("topright", designMatrix$Name, text.col = col, bty = "n", cex = 0.5)
  
  # Filter genes with low expression
  
  rawCount <- dim(countMatrix)
  
  keep.exprs <- filterByExpr(countMatrix, group = countMatrix$samples$group, lib.size = countMatrix$samples$lib.size)
  countMatrix <- countMatrix[keep.exprs,, keep.lib.sizes = FALSE] %>%
    calcNormFactors() 
  
  filterCount <- dim(countMatrix)
  
  glue::glue("{100 - round((filterCount[1]/rawCount[1])*100)}% of genes were filtered from {rawCount[2]} samples, \\
             where there were {rawCount[1]} genes before filtering and {filterCount[1]} genes after filtering")
  
  # Filtered density plot of log-CPM values 
  logCPM <- cpm(countMatrix, log = TRUE)
  plot(density(logCPM[,1]), col = col[1], lwd = 2, las =2 , main = "", xlab = "")
  title(main = "B. Filtered data", xlab = "Log-cpm")
  abline(v = logCPM.cutoff, lty = 3)
  for (i in 2:nsamples){
    den <- density(logCPM[,i])
    lines(den$x, den$y, col = col[i], lwd = 2)
  }
  legend("topright", designMatrix$Name, text.col = col, bty = "n", cex = 0.5)
  dev.off()
  
  # Interactive MDS plot
  Glimma::glMDSPlot(countMatrix,
                    groups = designMatrix,
                    path = getwd(),
                    folder = "interactiveMDS",
                    html = glue::glue("{tissueName}_MDS-Plot"),
                    launch = FALSE)
  
  # # MDS of all interactions
  # group <- interaction(designMatrix$Treatment, designMatrix$Sex)
  # plotMDS(countMatrix, col = as.numeric(group))
  
  # MDS of treatments simplified
  pdf(glue::glue("{tissueName}_treatmentMDS.pdf"), height = 8.5, width = 11)
  plotMDS(countMatrix,
          col = (as.numeric(countMatrix$samples$group)),
          main = glue::glue("{tissueName}_Treatment MDS"))
  dev.off()
  
  # Surrogate variables analysis --------------------------------------------
  
  # # Create model matrices, with null model for svaseq, and don't force a zero intercept
  # mm <- model.matrix(~Treatment + Sex + Litter,
  #                    data = designMatrix)
  # 
  # mm0 <- model.matrix(~1 + Sex + Litter,
  #                     data = designMatrix)
  # 
  # # svaseq requires normalized data that isn't log transformed
  # cpm <- cpm(countMatrix, log = FALSE)
  # 
  # # Calculate number of surrogate variables
  # nSv <- num.sv(cpm,
  #               mm,
  #               method = "leek")
  # 
  # # Estimate surrogate variables
  # svObj <- svaseq(cpm,
  #                 mm,
  #                 mm0,
  #                 n.sv = nSv)
  # 
  # # Update model to include surrogate variables
  # mm <- model.matrix(~Treatment + Sex + svObj$sv,
  #                    data = designMatrix)
  
  # Voom transformation and calculation of variance weights -----------------
  
  # Design
  mm <- model.matrix(~Treatment + Sex,
                     data = designMatrix)
  
  # Voom
  pdf(glue::glue("{tissueName}_voom_mean-variance_trend.pdf"), height = 8.5, width = 11)
  voomLogCPM <- voom(countMatrix,
                     mm,
                     plot = T)
  dev.off()
  
  # Make litter a random effect, since limma warns "coefficients not estimable" for some litters
  correlations <- duplicateCorrelation(voomLogCPM,
                                       mm,
                                       block = designMatrix$Litter)
  
  # Extract intraclass correlation within litters
  correlations <- correlations$consensus.correlation
  
  # Boxplots of logCPM values before and after normalization
  pdf(glue::glue("{tissueName}_normalization_boxplots.pdf"), height = 8.5, width = 11)
  par(mfrow=c(1,2))
  
  boxplot(logCPM, las = 2, col = col, main = "")
  title(main = "A. Unnormalised data", ylab = "Log-cpm")
  
  boxplot(voomLogCPM$E, las = 2, col = col, main = "")
  title(main = "B. Normalised data", ylab = "Log-cpm")
  
  dev.off()
  
  # Fitting linear models in limma ------------------------------------------
  
  # Wieght standard errors of log fold changes by within litter correlation 
  fit <- lmFit(voomLogCPM,
               mm,
               correlation = correlations,
               block = designMatrix$Litter)
  
  head(coef(fit))
  
  # Save normalized expression values for WGCNA
  voomLogCPM$E %>%
    write.xlsx(glue::glue("{tissueName}_voomLogCPMforWGCNA.xlsx"))
  
  
  # Split by contrast -------------------------------------------------------
  
  for (dose in 2:4){
    doseName = case_when(dose == 2 ~ 0.1,
                     dose == 3 ~ 1,
                     dose == 4 ~ 6)
    
    # Create DEG tibble -------------------------------------------------------
    
    print(glue::glue("Creating DEG list of {tissueName} samples for {doseName} PCB dosage"))
    
    DEGs <- fit %>%
      contrasts.fit(coef = dose) %>% # Change for different models (2,3,4)
      eBayes() %>%
      topTable(sort.by = "P", n = Inf) %>%
      rownames_to_column() %>% 
      tibble::as_tibble() %>%
      dplyr::rename(ensembl = rowname) %>% 
      dplyr::inner_join(grcm38, by = c("ensembl" = "ensgene")) %>% 
      dplyr::select(symbol, logFC, P.Value, adj.P.Val, ensembl, description) %>%
      dplyr::filter(P.Value < 0.05)
    
    # HTML report -------------------------------------------------------------
    
    print(glue::glue("Saving html report of {tissueName} samples for {doseName} PCB dosage"))
    
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
      write(glue::glue("{tissueName}_{doseName}_DEGs.html")) 
    
    # Plots -------------------------------------------------------------------
    
    # Heatmap
    
    print(glue::glue("Plotting heatmap of {tissueName} samples for {doseName} PCB dosage"))
    
    heatSamples <- designMatrix %>%
      dplyr::filter(Treatment == c("0", doseName)) %>%
      dplyr::select(Name) %>%
      purrr::flatten_chr()
    
    heatDesign <- designMatrix %>%
      dplyr::filter(Treatment == c("0", doseName))
    
    heatMatrix <- voomLogCPM$E[which(rownames(voomLogCPM$E) %in% DEGs$ensembl),] %>%
      tibble::as_tibble() %>%
      dplyr::select(one_of(heatSamples)) %>%
      as.matrix() # %>%
    #sweep(., 1, rowMeans(.))
    
    # https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
    gg_color_hue <- function(n = n){
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    
    gg_color <- c(gg_color_hue(length(levels(as.factor(heatDesign$Treatment)))))
    
    ColSideColors <- plyr::mapvalues(as.factor(heatDesign$Treatment),
                                     from = levels(as.factor(heatDesign$Treatment)),
                                     to = unique(gg_color)) 
    
    pdf(glue::glue("{tissueName}_{doseName}_heatmap.pdf"), height = 8.5, width = 11)
    
    heatmap.2(heatMatrix,
              scale = "row",
              labCol = heatDesign$Treatment,
              labRow = NA,
              col = rev(brewer.pal(11, name = "RdBu")),
              trace = "none",
              main = glue::glue("{nrow(DEGs)} Differentially Expressed Genes"),
              #key.xlab = "Z-score (log(cpm) - mean)",
              Rowv= as.dendrogram(hclust(dist(heatMatrix))),
              Colv = T,
              ColSideColors = as.character(ColSideColors))
    
    dev.off()
    
    # # MD plot
    # 
    # dt <- fit %>% decideTests() %>% summary()
    # 
    # pdf("MDplot.pdf", height = 8.5, width = 11)
    # 
    # plotMD(fit,
    #        column = 2,
    #        status = dt[,2],
    #        main = colnames(fit)[2])
    # 
    # dev.off()
    
    # Ontologies and Pathways -------------------------------------------------
    
    print(glue::glue("Performing GO and pathway analysis of {tissueName} samples for {doseName} PCB dosage"))
    
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
    
    write.xlsx(GO, file = glue::glue("{tissueName}_{doseName}_enrichr.xlsx"), sep = "")
    
    # Plot
    GOplot <- rbind(GO$GO_Biological_Process_2018[c(1:5),],
                    GO$GO_Cellular_Component_2018[c(1:5),],
                    GO$GO_Molecular_Function_2018[c(1:5),],
                    GO$KEGG_2016[c(1:5),]) %>%
      dplyr::as_tibble() %>%
      cbind(
        dplyr::as_tibble(
          c(
            rep("Biological Process", 5),
            rep("Cellular Component", 5),
            rep("Molecular Function", 5),
            rep("KEGG", 5)
          )
        )
      ) %>%
      dplyr::select(Term, P.value, value, Combined.Score) %>%
      dplyr::filter(P.value <= 0.05) %>%
      dplyr::mutate(P.value = -log10(P.value)) %>%
      dplyr::rename(`-log10.p-value` = P.value) %>%
      dplyr::rename(Database = value) %>% 
      dplyr::mutate(Database = stringr::str_replace(.$Database, "KEGG", "Pathway (KEGG)")) %>% 
      dplyr::mutate(Term = stringr::str_replace(.$Term, "\\(.*", "")) %>%
      dplyr::mutate(Term = stringr::str_replace(.$Term, "_.*", "")) %>%
      dplyr::mutate(Term = stringr::str_trim(.$Term)) %>%
      dplyr::mutate(Term = stringr::str_to_title(.$Term)) %>%
      dplyr::mutate(Database = factor(.$Database)) %>% 
      dplyr::mutate(Term = factor(.$Term, levels = unique(.$Term[order(forcats::fct_rev(.$Database), .$`-log10.p-value`)]))) %>% 
      ggplot2::ggplot(aes(x = Term, y = `-log10.p-value`, fill = Database, group = Database)) +
      geom_bar(stat = "identity", position = position_dodge(), color = "Black") +
      coord_flip() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            strip.text = element_text(size = 12)) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(y = expression("-log"[10](p))) +
      theme(axis.title.y = element_blank())
    
    ggsave(glue::glue("{tissueName}_{doseName}_enrichr_plot.pdf"),
           plot = GOplot,
           device = NULL,
           height = 8.5,
           width = 12)

  }
}
