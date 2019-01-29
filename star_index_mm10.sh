#!/bin/bash
#
#SBATCH --job-name=star_index
#SBATCH --workdir /share/lasallelab/Ben/PEBBLES/tag-seq
#SBATCH --ntasks=30 # Number of cores/threads
#SBATCH --mem=64000 # Ram in Mb
#SBATCH --partition=production 
#SBATCH --time=2-00:00:00

##########################################################################################
# Author: Ben Laufer
# Email: blaufer@ucdavis.edu 
# Modified from: https://leonjessen.wordpress.com/2014/12/01/how-do-i-create-star-indices-using-the-newest-grch38-version-of-the-human-genome/
##########################################################################################

###################
# Run Information #
###################

start=`date +%s`

hostname

THREADS=${SLURM_NTASKS}
MEM=$(expr ${SLURM_MEM_PER_CPU} / 1024)

echo "Allocated threads: " $THREADS
echo "Allocated memory: " $MEM

################
# Load Modules #
################

module load star/2.6.1d

###########################
# Download Genome (FASTA) #
###########################
# Must have un-placed/un-localized scaffolds, so use .dna.primary.assembly
# Section 2.2.1 https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

mkdir -p data/GRCm38/sequence

cd data/GRCm38/sequence/

wget ftp://ftp.ensembl.org/pub/release-95/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

cd ../../../

########################
# Download Genes (GTF) #
########################

mkdir -p data/GRCm38/annotation

cd data/GRCm38/annotation/

wget ftp://ftp.ensembl.org/pub/release-95/gtf/mus_musculus/Mus_musculus.GRCm38.95.gtf.gz

gunzip Mus_musculus.GRCm38.95.gtf.gz

cd ../../../

####################
# Build Star Index #
####################
# the splice-junction-data-base-overhang parameter should have a value of read length – 1

mkdir -p data/GRCm38/star_100/

call="STAR \
--runThreadN 30 \
--runMode genomeGenerate \
--genomeDir data/GRCm38/star_100/ \
--genomeFastaFiles data/GRCm38/sequence/Mus_musculus.GRCm38.dna.primary_assembly.fa \
--sjdbGTFfile data/GRCm38/annotation/Mus_musculus.GRCm38.95.gtf \
--sjdbOverhang 99"

echo $call
eval $call

###################
# Run Information #
###################

end=`date +%s`
runtime=$((end-start))
echo $runtime
