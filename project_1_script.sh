#!/bin/bash

echo "####################################################################################"
echo "####################################################################################"
echo “Running Zombie Tracker Pro  v1.0”
echo “Group 6”
echo "This program is for Research Purposes only"
echo "NOT for Diagnostic Use"
echo "####################################################################################"
echo "####################################################################################"
wait

# Create directories
mkdir ~/Project_1
mkdir ~/Project_1/RawData
mkdir ~/Project_1/ProcessedData
mkdir ~/Project_1/ProcessedData/BAM
mkdir ~/Project_1/ProcessedData/VCF
mkdir ~/Project_1/ProcessedData/MFA
mkdir ~/Project_1/ProcessedData/FASTA

# Copy reference genome FASTA to RawData folder
cp /projects/micb405/resources/project_1/ref_genome.fasta ~/Project_1/RawData

# Index reference genome and move reference files to ref_index folder
bwa index ~/Project_1/RawData/ref_genome.fasta

echo "Running BWA analysis for fastq files"
wait

# ======================================================
# ======== Alignment Steps and Variant Calling =========
# ======================================================
for fastq in /projects/micb405/resources/project_1/*_1.fastq.gz \
; do prefix=$(basename $fastq | sed 's/_1.fastq.gz//g') \
; echo "Analyzing $prefix"; bwa mem -t 20 ~/Project_1/RawData/ref_genome.fasta /projects/micb405/resources/project_1/$prefix\_1.fastq.gz /projects/micb405/resources/project_1/$prefix\_2.fastq.gz | samtools view -b - >~/Project_1/ProcessedData/BAM/$prefix.bam \
; echo "Sorting $prefix"; samtools sort -@ 20 ~/Project_1/ProcessedData/BAM/$prefix.bam -o ~/Project_1/ProcessedData/BAM/$prefix.sorted.bam \
; samtools rmdup ~/Project_1/ProcessedData/BAM/$prefix.sorted.bam ~/Project_1/ProcessedData/BAM/$prefix.sorted.rmdup.bam \
; echo "Indexing $prefix"; samtools index ~/Project_1/ProcessedData/BAM/$prefix.sorted.rmdup.bam \
; echo "Calling Variants for $prefix"; bcftools mpileup --fasta-ref ~/Project_1/RawData/ref_genome.fasta ~/Project_1/ProcessedData/BAM/$prefix.sorted.rmdup.bam | bcftools call -mv - >~/Project_1/ProcessedData/VCF/$prefix.raw.vcf \
; echo "Filtering Variants for $prefix"; bcftools filter --threads 12 --exclude "QUAL < 200" ~/Project_1/ProcessedData/VCF/$prefix.raw.vcf; done
	

# ==================================
# == Multiple Sequence Alignment ===
# ==================================
echo "####################################################################################"
echo "Converting VCF files into FASTA file called snps_200r.fasta"
echo "####################################################################################"
python /projects/micb405/resources/vcf_to_fasta_het.py -x ~/Project_1/ProcessedData/VCF/ snps_200r
cp ~/Project_1/ProcessedData/VCF/snps_200r* ~/Project_1/ProcessedData/FASTA/

echo "####################################################################################"
echo "Running Muscle for FASTA files"
echo "####################################################################################"
muscle -in ~/Project_1/ProcessedData/FASTA/snps_200r.fasta -out ~/Project_1/ProcessedData/MFA/snps_200r_muscle.mfa

echo "Running Trimal for MUSCLE files"
trimal -automated1 -in ~/Project_1/ProcessedData/MFA/snps_200r_muscle.mfa -out ~/Project_1/ProcessedData/MFA/snps_200r_muscle.trimal.mfa

# ==================================
# == Building phylogenetic trees ===
# ==================================
echo "####################################################################################"
echo "Running RAxML for MFA files"
echo "####################################################################################"
raxml-ng --all --msa ~/Project_1/ProcessedData/MFA/snps_200r_muscle.trimal.mfa --model GTR+G4 --tree rand{100} --bs-trees 20 --threads 2 --seed 12345

echo "####################################################################################"
echo "Running FastTree for MFA files"
echo "####################################################################################"
FastTree ~/Project_1/ProcessedData/MFA/snps_200r_muscle.trimal.mfa 1>~/Project_1/ProcessedData/MFA/snps_200r_muscle.fastTree.nwk

echo "####################################################################################"
echo "####################################################################################"
echo "Finished Running Tree Generator"
echo “Zombie Transmissions can now be identified” 
echo "####################################################################################"
echo "####################################################################################"
