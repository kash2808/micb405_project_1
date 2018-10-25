#!/bin/bash

# Create Project_1 folder, copy ref_genome.fasta to RawData folder
mkdir Project_1
mkdir ~/Project_1/RawData
cp /projects/micb405/resources/Project_1/ref_genome.fasta ~/Project_1/RawData

# Create ProcessedData, BWA_output, and ref_index folder
mkdir ~/Project_1/ProcessedData
mkdir ~/Project_1/ProcessedData/BWA_output
mkdir ~/Project_1/ProcessedData/BWA_output/ref_index

# Index reference genome and move reference files to ref_index folder
echo "Generating index"
bwa index ~/Project_1/RawData/ref_genome.fasta
mv ref* ~/Project_1/ProcessedData/BWA_output/ref_index

# ==================================
# ======== Alignment Steps =========
# ==================================
# Align each sequence to the reference genome index and create .bam files

for fastq in /projects/micb405/resources/Project_1/*_1.fastq.gz;
do
  prefix=$(basename $fastq | sed 's/_1.fastq.gz//g')
  echo "Analyzing $prefix"
  bwa mem -t 10 ~/Project_1/RawData/ref_genome.fasta \
  /projects/micb405/resources/Project_1/$prefix\_1.fastq.gz /projects/micb405/resources/Project_1/$prefix\_2.fastq.gz | \
  samtools view -b - >~/Project_1/ProcessedData/BAM/$prefix.bam
done


# ==================================
# ======= Variant Calling ==========
# ==================================
mkdir ~/Project_1/ProcessedData/BAM
mkdir ~/Project_1/ProcessedData/VCF

for bam in ~/Project_1/ProcessedData/BAM/*bam;
do
	prefix=$(basename $bam | ’s/.bam//g’)
	echo "Analyzing $prefix"
	# Sorts bam file by alignment position
	samtools sort ~/Project_1/ProcessedData/BAM/$prefix.bam \
	1>~/Project_1/ProcessedData/BAM/$prefix.sorted.bam

	# Remove duplicates
	samtools rmdup \
	~/Project_1/ProcessedData/BAM/$prefix.sorted.bam \
	~/Project_1/ProcessedData/BAM/$prefix.sorted.rmdup.bam

	# Index bam file
	samtools index ~/Project_1/ProcessedData/BAM/$prefix.sorted.rmdup.bam

	# Variant calling
	bcftools mpileup --fasta-ref ~/Project_1/RawData/ref_genome.fasta /
	~/Project_1/ProcessedData/BAM/$prefix.sorted.rmdup.bam | bcftools call -mv - \
	>~/Project_1/ProcessedData/VCF/$prefix.VS.raw.vcf

	# Variant Filtering
	bcftools filter --exclude "QUAL < 200" /Project_1/ProcessedData/VCF/$prefix.VS.raw.vcf
done 


# ==================================
# == Multiple Sequence Alignment ===
# ==================================

mkdir ~/Project_1/ProcessedData/MFA

# Dr. Jennifer Gardy’s script to convert .vcf to .fasta files
python /projects/micb405/resources/vcf_to_fasta_het.py -x ~/Project_1/ProcessedData/VCF/spooky.snps

# MSA with MUSCLE
muscle -in ~/Project_1/ProcessedData/VCF/spooky.snps.fasta -out ~/Project_1/ProcessedData/MFA/spooky.snps_muscle.mfa
 
# MSA with MAFFT
mafft ~/Project_1/ProcessedData/VCF/spooky.snps.fasta > ~/Project_1/ProcessedData/MFA/spooky.snps_mafft.mfa

# Trimal for MUSCLE file
trimal -automated1 \
-in ~/Project_1/ProcessedData/MFA/spooky.snps_muscle.mfa \
-out ~/Project_1/ProcessedData/MFA/spooky.snps_muscle.trimal.mfa

# Trimal for MAFFT File
trimal -automated1 \
-in ~/Project_1/ProcessedData/MFA/spooky.snps_mafft.mfa \
-out ~/Project_1/ProcessedData/MFA/spooky.snps_mafft.trimal.mfa


# ==================================
# == Building phylogenetic trees ===
# ==================================

# RAxML for maximum-likelihood for MUSCLE
raxml-ng --all --msa ~/Project_1/ProcessedData/MFA/spooky.snps_muscle.mfa \
--model LG+G2 --tree rand{100} --bs-trees 20 \
--threads 12 --seed 12345

# FastTree for MUSCLE
FastTree ~/Project_1/ProcessedData/MFA/spooky.snps_muscle.mfa \
1>~/Project_1/ProcessedData/MFA/spooky.snps_muscle.nwk

# RAxML for maximum-likelihood for MAFFT
raxml-ng --all --msa ~/Project_1/ProcessedData/MFA/spooky.snps_mafft.mfa \
--model LG+G2 --tree rand{100} --bs-trees 20 \
--threads 12 --seed 12345

# FastTree for MAFFT
FastTree ~/Project_1/ProcessedData/MFA/spooky.snps_mafft.mfa \
1>~/Project_1/ProcessedData/MFA/spooky.snps_mafft.nwk