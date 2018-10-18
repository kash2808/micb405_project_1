#!/bin/bash

# Create project_1, RawData folders and copy ref_genome.fasta and all fastq.gz files
mkdir project_1
mkdir project_1/RawData
cp /projects/micb405/resources/project_1/ref_genome.fasta ~/project_1/RawData
cp /projects/micb405/resources/project_1/*.fastq.gz ~/project_1/RawData
cd project_1/RawData

# Unzip fastq files
gunzip *.fastq.gz

# Create ProcessedData, BWA_output, and ref_index folder
mkdir ~/project_1/ProcessedData
mkdir ~/project_1/ProcessedData/BWA_output
mkdir ~/project_1/ProcessedData/BWA_output/ref_index

# Index reference genome and move reference files to ref_index folder
bwa index ref_genome.fasta
mv ref* ~/project_1/ProcessedData/BWA_output/ref_index

# ========================
# === Alignment Steps ====
# ========================
clear

echo "Running BWA mem"

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Patient_1_1.fastq ~/project_1/RawData/Patient_1_2.fastq 1>~/project_1/ProcessedData/BWA_output/Patient_1.sam 2>~/project_1/ProcessedData/BWA_output/Patient_1_log.txt 

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Patient_2_1.fastq ~/project_1/RawData/Patient_2_2.fastq 1>~/project_1/ProcessedData/BWA_output/Patient_2.sam 2>~/project_1/ProcessedData/BWA_output/Patient_2_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Patient_3_1.fastq ~/project_1/RawData/Patient_3_2.fastq 1>~/project_1/ProcessedData/BWA_output/Patient_3.sam 2>~/project_1/ProcessedData/BWA_output/Patient_3_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Patient_4_1.fastq ~/project_1/RawData/Patient_4_2.fastq 1>~/project_1/ProcessedData/BWA_output/Patient_4.sam 2>~/project_1/ProcessedData/BWA_output/Patient_4_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Patient_5_1.fastq ~/project_1/RawData/Patient_5_2.fastq 1>~/project_1/ProcessedData/BWA_output/Patient_5.sam 2>~/project_1/ProcessedData/BWA_output/Patient_5_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Patient_6_1.fastq ~/project_1/RawData/Patient_6_2.fastq 1>~/project_1/ProcessedData/BWA_output/Patient_6.sam 2>~/project_1/ProcessedData/BWA_output/Patient_6_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Patient_7_1.fastq ~/project_1/RawData/Patient_7_2.fastq 1>~/project_1/ProcessedData/BWA_output/Patient_7.sam 2>~/project_1/ProcessedData/BWA_output/Patient_7_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Patient_8_1.fastq ~/project_1/RawData/Patient_8_2.fastq 1>~/project_1/ProcessedData/BWA_output/Patient_8.sam 2>~/project_1/ProcessedData/BWA_output/Patient_8_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Patient_9_1.fastq ~/project_1/RawData/Patient_9_2.fastq 1>~/project_1/ProcessedData/BWA_output/Patient_9.sam 2>~/project_1/ProcessedData/BWA_output/Patient_9_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Patient_10_1.fastq ~/project_1/RawData/Patient_10_2.fastq 1>~/project_1/ProcessedData/BWA_output/Patient_10.sam 2>~/project_1/ProcessedData/BWA_output/Patient_10_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Patient_11_1.fastq ~/project_1/RawData/Patient_11_2.fastq 1>~/project_1/ProcessedData/BWA_output/Patient_11.sam 2>~/project_1/ProcessedData/BWA_output/Patient_11_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Patient_12_1.fastq ~/project_1/RawData/Patient_12_2.fastq 1>~/project_1/ProcessedData/BWA_output/Patient_12.sam 2>~/project_1/ProcessedData/BWA_output/Patient_12_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Patient_13_1.fastq ~/project_1/RawData/Patient_13_2.fastq 1>~/project_1/ProcessedData/BWA_output/Patient_13.sam 2>~/project_1/ProcessedData/BWA_output/Patient_13_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Patient_14_1.fastq ~/project_1/RawData/Patient_14_2.fastq 1>~/project_1/ProcessedData/BWA_output/Patient_14.sam 2>~/project_1/ProcessedData/BWA_output/Patient_14_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Patient_15_1.fastq ~/project_1/RawData/Patient_15_2.fastq 1>~/project_1/ProcessedData/BWA_output/Patient_15.sam 2>~/project_1/ProcessedData/BWA_output/Patient_15_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Patient_16_1.fastq ~/project_1/RawData/Patient_16_2.fastq 1>~/project_1/ProcessedData/BWA_output/Patient_16.sam 2>~/project_1/ProcessedData/BWA_output/Patient_16_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Bat_1.fastq ~/project_1/RawData/Bat_2.fastq 1>~/project_1/ProcessedData/BWA_output/Bat.sam 2>~/project_1/ProcessedData/BWA_output/Bat_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Cat_1.fastq ~/project_1/RawData/Cat_2.fastq 1>~/project_1/ProcessedData/BWA_output/Cat.sam 2>~/project_1/ProcessedData/BWA_output/Cat_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Guinea_pig_1.fastq ~/project_1/RawData/Guinea_pig_2.fastq 1>~/project_1/ProcessedData/BWA_output/Guinea_pig.sam 2>~/project_1/ProcessedData/BWA_output/Guinea_pig_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Rabid_raccoon_1_1.fastq ~/project_1/RawData/Rabid_raccoon_1_2.fastq 1>~/project_1/ProcessedData/BWA_output/Rabid_raccoon_1.sam 2>~/project_1/ProcessedData/BWA_output/Rabid_raccoon_1_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Rabid_raccoon_2_1.fastq ~/project_1/RawData/Rabid_raccoon_2_2.fastq 1>~/project_1/ProcessedData/BWA_output/Rabid_raccoon_2.sam 2>~/project_1/ProcessedData/BWA_output/Rabid_raccoon_2_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Rabid_raccoon_3_1.fastq ~/project_1/RawData/Rabid_raccoon_3_2.fastq 1>~/project_1/ProcessedData/BWA_output/Rabid_raccoon_3.sam 2>~/project_1/ProcessedData/BWA_output/Rabid_raccoon_3_log.txt

bwa mem -10 ~/project_1/ProcessedData/BWA_output/ref_index ~/project_1/RawData/Rabid_raccoon_4_1.fastq ~/project_1/RawData/Rabid_raccoon_4_2.fastq 1>~/project_1/ProcessedData/BWA_output/Rabid_raccoon_4.sam 2>~/project_1/ProcessedData/BWA_output/Rabid_raccoon_4_log.txt

echo "BWA mem completed"
