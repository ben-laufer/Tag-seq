# Tag-seq
### A 3’ Tag RNA-seq repository for the alignment and statistical analysis of gene expression data

The Tag-seq repository consists of two main pipelines, which are not yet finished. The first pipeline, which utilizes Unix shell scripts, takes you from raw fastq files to gene counts matrices for samples prepared using the QuantSeq 3' mRNA-Seq Library Prep Kit FWD for Illumina (Lexogen) and sequenced with single end reads. These shell scripts handle the indexing of the reference genome and gene annotations, trimming, alignment, extraction of gene counts, and generate a quality control report. The second pipeline, which is an R script, takes the gene count matrices and performs preprocessing, normalization, differential gene expression (DGE) analyses, annotations, gene ontology and enrichment testing, and data visualization.

**Note:** This pipeline is based on old data without Unique Molecular Identifiers (UMIs), which are now becoming a standard for newer datasets. Please check if your data has them, since they require an extra step of processing that is not included in this repository. 

### Table of Contents

1. [Installation](https://github.com/ben-laufer/Tag-seq#installation)
2. [Indexing](https://github.com/ben-laufer/Tag-seq#indexing)
3. [Single End (SE) Sequencing](https://github.com/ben-laufer/Tag-seq#single-end-se-sequencing)
   1. [Trimming](https://github.com/ben-laufer/Tag-seq#trimming)
   2. [Alignment](https://github.com/ben-laufer/Tag-seq#alignment)
4. [Quality Control](https://github.com/ben-laufer/Tag-seq#quality-control)
4. [Statistical Analysis](https://github.com/ben-laufer/Tag-seq#statistical-analysis)
6. [Citation](https://github.com/ben-laufer/Tag-seq#citation)
7. [Acknowledgements](https://github.com/ben-laufer/Tag-seq#acknowledgements)

## Installation

This workflow utilizes the following packages, which need to be installed and in your path:
1. [FastQC](https://github.com/s-andrews/FastQC)
2. [BBMap](sourceforge.net/projects/bbmap/)
3. [STAR](https://github.com/alexdobin/STAR)
4. [Samtools](http://www.htslib.org)
5. [MultiQC](http://multiqc.info)

I recommend using [Bioconda](https://bioconda.github.io) to install and manage the package updates, which can be accomplished by:

`conda install -c bioconda fastqc bbmap star samtools multiqc`

## Indexing

Both the genome (FASTA) and gene annotations (GTF) are required to create the STAR indexes. An example of how to download the required files and create the index is provided in the `star_index_mm10.sh` script. The indexes are specific to both species and read length.

The complete genome folder structure should appear as:

```
├── data
│   ├── GRCm38
│   │   ├── star_100
│   │   ├── sequence
│   │   ├── annotation
│   ├── truseq_rna.fa.gz
```
For convenience the truseq_rna.fa.gz from BBMap is provided in this repository. 

## Single End (SE) Sequencing

The primary workflow is `tag-seq.sh`, which is a modified version of the [Lexogen example](https://www.lexogen.com/quantseq-data-analysis/). The modifications are two-fold:
1.  The script is optimized to run on a HPCC (High-Performance Computing Cluster) with the SLURM workload manager and offers more automation of file handling.
2. Additional files to aid with quality control and statistical analysis are also generated.

This pipeline requires a specific directory structure:

1.	Create a parent directory for the project
2.	Within that parent project directory, add a text file called “task_samples.txt”, where each new line contains the entire sample name exactly as it appears on the fastq read pair files, aside from the file extension (“.fastq.gz”). Also, if you’re using excel or a windows desktop, you will need to change the linebreaks from windows to unix, which can be done using BBEdit or:
`awk '{ sub("\r$", ""); print }' task_samples_windows.txt > task_samples.txt`
3.	Within that parent directory create a folder called “raw_sequences” that contains all the raw fastq files (.fastq.gz)

Overall, the directory tree structure should be the following:

```
├── Project
│   ├── raw_sequences
│   │   ├── sample1.fastq.gz
│   │   ├── sample2.fastq.gz
│   ├── task_samples.txt
│   ├── data
```

The script can be executed using the following command from the working directory:

`sbatch --array=1-87 tag-seq.sh`

The `--array` argument represents which samples in the list to run. Here we are running samples 1 to 87. You could run select samples using the following format --array=2,4-12

### Trimming

Trimming is performed to remove adapter contamination, polyA tail read through, and low quality tails. The random priming approach also biases the first 12 bp of sequence; however, this does not need to be manually trimmed if using a local aligner like STAR (read more [here](https://dnatech.genomecenter.ucdavis.edu/tag-seq-gene-expression-profiling/)).

### Alignment

After trimming, alignment and the generation of gene counts are performed in a single step that is then followed by indexing of the BAM files.

## Quality Control

The `tag-seq-QC.sh` script is the final step in the process, which a generates a html report and places the gene count matrices into a single folder.

## Statistical Analysis

The `tag-seq.R` script is intended for the statistical analysis of the gene count matrices with the goal of identifying differentially expressed genes. It utilizes Limma-Voom and integrates two tutorials (see [1](https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html) and [2](https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html)). The script also offers additional custom functions related to statistical analysis, annotation, and data visualization. 

In this linear model, sex is modeled as a fixed effect. Litter is modeled as a random effect for two reasons. The first reason is because limma warns "coefficients not estimable" for some litters when modeling them as a fixed effect. The second reason is due to the [nested design](https://support.bioconductor.org/p/11956/), where the treatment is applied to the dam but the effects are measured in multiple offspring from her litter. Additionally, surrogate variable analysis is not utilized, as it has been shown to exaggerate group differences in similar designs, where it was shown that using a blocking factor is a more statistically appropriate approach (read more in [1](https://www.ncbi.nlm.nih.gov/pubmed/26272994) and [2](https://www.ncbi.nlm.nih.gov/pubmed/27780809)).

There is also a rough version of a weighted gene co-expression network analysis (WGCNA) script that I've modified from the [tutorials](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html) and used for a few DNA methylation and gene expression projects, however, it isn't yet complete in terms of automation, error checking, and minor bugs. Furthermore, it appears that WGCNA is not appropriate for modelling non-monotonic responses and there are [not enough samples](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html) to analyze each dose individually in this project.

## Citation

If you use **Tag-seq** in published research please cite this repository.

## Acknowledgements

 I would like to thank Blythe Durbin-Johnson for statistical consulting and Annie Vogel Ciernia for invaluable discussions related to the bioinformatic and statistical approaches utilized in this repository. I would also like to thank [Matt Settles](https://github.com/msettles) from the [UC Davis Bioinformatics Core](https://github.com/ucdavis-bioinformatics) for [examples of tidy code](https://github.com/ucdavis-bioinformatics-training/A-Primer-on-Using-the-Bioinformatics-Core-Administrated-Servers-and-Cluster-s-/tree/master/examples) and introducing me to the tidyverse. Finally, I would like to thank [Ian Korf](https://github.com/KorfLab) for discussions related to some of the bioinformatic approaches utilized in this repository. 
