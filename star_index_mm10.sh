#!/bin/bash
#
#SBATCH --job-name=star_index
#SBATCH --workdir /share/lasallelab/Ben/PEBBLES/tag-seq
#SBATCH --ntasks=20 # Number of cores/threads
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

mkdir -p data/GRCm38/sequence

cd data/GRCm38/sequence/

wget ftp://ftp.ensembl.org/pub/release-95/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.{1..19}.fa.gz

wget ftp://ftp.ensembl.org/pub/release-95/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.{MT,X,Y}.fa.gz

gunzip -c Mus_musculus.GRCh38.dna.chromosome.* > GRCm38_r95.all.fa

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

mkdir -p data/GRCm38/star_indices_overhang100

# the splice-junction-data-base-overhang parameter should have a value of read length – 1
call="star \
--runThreadN 20 \
--runMode genomeGenerate \
--genomeDir data/GRCh38/star_indices_overhang100/ \
--genomeFastaFiles data/GRCh38/sequence/GRCm38_r95.all.fa \
--sjdbGTFfile data/GRCh38/annotation/Mus_musculus.GRCm38.95.gtf.gz \
--sjdbOverhang 99"

echo $call
eval $call

###################
# Run Information #
###################

end=`date +%s`
runtime=$((end-start))
echo $runtime
