# Tag-seq
### A 3’ Tag RNA-seq pipeline for the alignment of gene expression data

The Tag-seq pipeline takes you from raw fastq files to gene counts for samples prepared using the QuantSeq 3' mRNA-Seq Library Prep Kit FWD for Illumina (Lexogen) and sequenced with single end reads. The scripts cover the indexing of the reference genome and gene annotations, filtering, alignment, extraction of gene counts, quality control. The gene count matrices can then be utilized for statistical analysis in R. 

### Table of Contents

1. [Installation](https://github.com/ben-laufer/Tag-seq#installation)
2. [Indexing](https://github.com/ben-laufer/Tag-seq#indexing)
3. [Single End (SE) Sequencing](https://github.com/ben-laufer/Tag-seq#single-end-se-sequencing)
   1. [Filtering](https://github.com/ben-laufer/CpG_Me#filtering)
   2. [Alignment](https://github.com/ben-laufer/CpG_Me#alignment)
4. [Quality Control](https://github.com/ben-laufer/Tag-seq#qaulity-control)
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

### Filtering

Filtering is performed to remove adapter contamination, polyA tail read through, and low quality tails. The random priming approach also biases the first 12 bp of sequence; however, this does not need to be manually trimmed if using a local aligner like STAR (read more [here](https://dnatech.genomecenter.ucdavis.edu/tag-seq-gene-expression-profiling/)).

### Alignment

After filtering, alignment and the generation of gene counts are performed in a single step that is then followed by indexing of the BAM files.

## Quality Control

The `tag-seq-QC.sh` script is the final step in the process, which a generates a html report and places the gene count matrices into a single folder.

## Statistical Analysis

While the above scripts are stable, the `tag-seq.R` script is a work in progress intended for the statistical analysis of the gene count matrices. It utilizes [Limma-Voom](https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html) and [biobroom](https://www.bioconductor.org/packages/release/bioc/vignettes/biobroom/inst/doc/biobroom_vignette.html) for the analysis of differential gene expression.

## Citation

If you use **Tag-seq** in published research please cite this repository.

## Acknowledgements
I would like to thank [Matt Settles](https://github.com/msettles) from the [UC Davis Bioinformatics Core](https://github.com/ucdavis-bioinformatics) for [examples of tidy code](https://github.com/ucdavis-bioinformatics-training/A-Primer-on-Using-the-Bioinformatics-Core-Administrated-Servers-and-Cluster-s-/tree/master/examples). I would also like to thank Blythe Durbin-Johnson for statistical consulting and Annie Vogel Ciernia for invaluable discussions related to the bioinformatic and statistical approaches utilized in this repository. Finally, I would like to thank [Ian Korf](https://github.com/KorfLab) for invaluable discussions related to the bioinformatic approaches utilized in this repository. 

