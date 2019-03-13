# This file contains the commands used to rename and merge fastq files across multiple runs/lanes

# Test renaming of samples with duplicate names from different runs
rename -n 's/L001/L002/' *.fastq.gz 

# Use rename command by removing -n
rename 's/L001/L002/' *.fastq.gz

# Create file of unique IDs based on _ delimiter and first 3 strings
ls -1 *R1*.gz | awk -F '_' '{print $1"_"$2"_"$3}' | sort | uniq > ID

# Test merge command
for i in `cat ./ID`; do echo cat $i\_*.fastq.gz > $i\_merged_R1.fastq.gz; done

# Use merge command by removing echo
for i in `cat ./ID`; do cat $i\_*.fastq.gz > $i\_merged_R1.fastq.gz; done

# ref 1: https://stackoverflow.com/questions/1086502/rename-multiple-files-in-unix
# ref 2: https://www.biostars.org/p/317385/
