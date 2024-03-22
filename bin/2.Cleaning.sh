#!bin/bash

##Trimming of reads

#Create a folder to store trimmed data
mkdir -p ../data/trimmed

#Set the input and output paths
in_path="../data"
output_un_path="../data/unpaired"

#Filter the reads with a quality equal or above 30
#Create arrays with the name of the input (forward and reverse) files
for i in {1..18}; do
printf "\n" #Insert a line break
echo "t0$i" #Print the name of each sample
r1="t0$i-R1.fastq.gz"
r2="t0$i-R2.fastq.gz"
#Create arrays with the name of the output (paired and unpaired files)
r1p="t0$i-R1.paired.fastq.gz"
r2p="t0$i-R2.paired.fastq.gz"
r1u="t0$i-R1.unpaired.fastq.gz"
r2u="t0$i-R2.unpaired.fastq.gz"
#Run trimmomatic with AVGQUAL
trimmomatic PE $in_path/$r1 $in_path/$r2 $in_path/$r1p $output_un_path/$r1u $in_path/$r2p $output_un_path/$r2u AVGQUAL:25
done &> $in_path/read_trimming.txt

#Trim the Illumina adapters from the quality-filtered reads
#Generate a list with the name of the input (forward and reverse) files
for i in {1..18}; do
printf "\n"
echo "t0$i"
R1="t0$i-R1.paired.fastq.gz"
R2="t0$i-R2.paired.fastq.gz"
#List with the name of the output files (forward and reverse)
R1c="t0$i-R1.clean.paired.fastq.gz"
R2c="t0$i-R2.clean.paired.fastq.gz"
#Run cutadapt
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o $in_path/$R1c -p $in_path/$R2c $in_path/$R1 $in_path/$R2
done &> $in_path/adapter_trimming.txt

#Perform MINLENGTH of *.clean.fastq.gz reads
for i in {1..18}; do
printf "\n"
R1a="t0$i-R1.clean.fastq.gz"
R2a="t0$i-R2.clean.fastq.gz"
#Array for the output files
R1pa="t0$i-R1.clean.paired.fastq.gz"
R2pa="t0$i-R2.clean.paired.fastq.gz"
R1ua="t0$i-R1.clean.unpaired.fastq.gz"
R2ua="t0$i-R2.clean.unpaired.fastq.gz"
echo "#######SAMPLE"
echo "t0$i"
#Run trimmomatic with MINLEN
trimmomatic PE $out_p_path/$R1a $out_p_path/$R2a $out_p_path/$R1pa $output_un_path/$R1ua $out_p_path/$R2pa $output_un_path/$R2ua MINLEN:70
done &> $out_p_path/minlen4.txt


