#!/bin/bash
#
#SBATCH --job-name=tagQC
#SBATCH --workdir /share/lasallelab/Ben/PEBBLES/tag-seq
#SBATCH --ntasks=1 # Number of cores/threads
#SBATCH --mem=2000 # Ram in Mb
#SBATCH --partition=production 
#SBATCH --time=0-00:30:00

##########################################################################################
# Author: Ben Laufer
# Email: blaufer@ucdavis.edu 
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

module load multiqc/1.6

###########
# MultiQC #
###########

call="multiqc
. \
 --config multiqc_config.yaml"

echo $call
eval $call

########
# Copy #
########

mkdir bam_sorted
"$(find `.` -name '*Aligned.sortedByCoord.out.bam' -print0 | xargs -0 cp -t bam_sorted)"

###################
# Run Information #
###################

end=`date +%s`
runtime=$((end-start))
echo $runtime
