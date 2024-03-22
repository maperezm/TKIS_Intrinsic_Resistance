#!bin/bash

###RSEM Counts estimation 
#Create a directory to store the  Results 
mkdir -p ../rsem

for i in {1..18}; do
printf "\n"
echo "t0$i"?
#Generate the list with the name of the forward (r1) and reverse (r2) files
r1="t0$i-R1.clean.paired.fastq.gz"
r2="t0$i-R2.clean.paired.fastq.gz"
#execute RSEM Program
rsem-calculate-expression --star-gzipped-read-file \
    --star --paired-end -p 32 \
    ../data$r1.fastq.gz ../data$r2.fastq.gz
    ../index / \ 
    ../data/t0$i
done 