#!bin/bash

###FastQC analysis and percentage of Q30 reads 

#Create a directory "fastqc_output" inside the directory "OSF_Storage"
mkdir -p ../data/fastqc_output
mkdir -p ../multiquc_output
#Start a FastQC analysis of each sample and send the results (.html and .zip) to fastqc_output
for i in ../data/*.fastq.gz; do
fastqc -o ../data/fastqc_output $i
done

##Run multiqc
multiqc ../fastqc_output/*.zip -o ../multiquc_output
