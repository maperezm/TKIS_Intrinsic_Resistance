#!bin/bash

###Alignment of reads
#Create a directory to store the index 
mkdir -p ../index

#Create the index for the human genome
STAR --runThreadN 24 \ #Number of cores
    --runMode genomeGenerate \ #Option to generate the index
     --genomeDir ../index \ #Path to store the index
     --genomeFastaFiles ../genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \ #Path to the fasta file of the genome
     --sjdbGTFfile ../annotation/Homo_sapiens.GRCh38.83.gtf \ #Path to the gtf file of the annotation
    --sjdbOverhang 100

#Align the reads using STAR
for i in {1..18}; do # make aligment of all samples
printf "\n"
echo "t0$i"?
#Generate the list with the name of the forward (r1) and reverse (r2) files
r1="t0$i-R1.clean.paired.fastq.gz" #names asigment after trimming 
r2="t0$i-R2.clean.paired.fastq.gz"
STAR --genomeDir ../index \ #Path to the index previously generated
     --runThreadN 32 \ #Number of cores
     --readFilesIn ../data/$r1 ../data/$r2 \ #Path to the forward and reverse reads
     --outFileNamePrefix ../data/t0$i \ #Prefix for the output
     --outSAMtype BAM SortedByCoordinate \ #Specify the type of file in the output (BAM format)
     --quantMode TranscriptomeSAM
done 


