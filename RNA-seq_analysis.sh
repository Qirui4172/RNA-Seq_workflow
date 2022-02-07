#!/bin/bash

# Qirui Zhang (qirui.zhang@med.uni-greifswald.de)
# 01-02-2022


#=====================================================================================================
WorkDir="/home/nhuber/qirui/20210929_Combi-PR-MI503-Nuria/Analysis_qirui/"
cd ${WorkDir}
AdapterDir="/home/nhuber/Apps/Trimmomatic-0.39/adapters"
GenomeIndx="/home/nhuber/Genomes/hg38/STAR/"
Gtf="/home/nhuber/Genomes/hg38/hg38.ncbiRefSeq.gtf"
Samples=("Com1" "Com2" "Com3" "Com4" "D1" "D2" "D3" "D4" "MI1" "MI2" "MI3" "MI4")


#=====================================================================================================
# Step0: Create folders

mkdir -p QC QC/FastqRaw QC/FastqTrimmed QC/AlignmentRaw QC/AlignmentFiltered QC/multiqc FastqRaw AlignmentRaw AlignmentFiltered ReadsCount


#=====================================================================================================
# Step1: Raw fastq quality control

for sample in ${Samples[@]};do
  time=$(date "+%Y-%m-%d %H:%M:%S")
  echo -e "$time quality control on raw fastq of $sample ..."

  fastqc -t 12 -o QC/FastqRaw/ -q FastqRaw/${sample}_*_R1_001.fastq.gz FastqRaw/${sample}_*_R2_001.fastq.gz
done


#=====================================================================================================
# Step2: Trim adapters

for sample in ${Samples[@]};do
  time=$(date "+%Y-%m-%d %H:%M:%S")
  echo -e "\n$time Trimming adapters of $sample ..."

  java -jar ~/Apps/Trimmomatic-0.39/trimmomatic-0.39.jar PE FastqRaw/${sample}_*_R1_001.fastq.gz FastqRaw/${sample}_*_R2_001.fastq.gz FastqTrimmed/${sample}_R1_paired.fastq.gz FastqTrimmed/${sample}_R1_unpaired.fastq.gz FastqTrimmed/${sample}_R2_paired.fastq.gz FastqTrimmed/${sample}_R2_unpaired.fastq.gz ILLUMINACLIP:${AdapterDir}/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -threads 12
done

# quality control
for sample in ${Samples[@]};do
  time=$(date "+%Y-%m-%d %H:%M:%S")
  echo -e "\n$time quality control of $sample ..."

  fastqc -t 12 -o QC/FastqTrimmed/ -q FastqTrimmed/${sample}_R1_paired.fastq.gz FastqTrimmed/${sample}_R2_paired.fastq.gz
done


#=====================================================================================================
# Step3: Mapping

for sample in ${Samples[@]};do
  time=$(date "+%H:%M:%S %Y-%m-%d")
  echo -e "\n$time Mapping $sample ..."

  # mapping
  mkdir -p AlignmentRaw/${sample}
  STAR --runThreadN 12 --genomeDir $GenomeIndx --readFilesIn FastqTrimmed/${sample}_R1_paired.fastq.gz FastqTrimmed/${sample}_R2_paired.fastq.gz --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix AlignmentRaw/${sample}/ --readFilesCommand zcat

  # quality control
  samtools flagstat -@ 12 AlignmentRaw/${sample}/Aligned.sortedByCoord.out.bam > QC/AlignmentRaw/${sample}_flagstat.txt
done


#=====================================================================================================
# Step4: Filter alignment

for sample in ${Samples[@]};do
  time=$(date "+%H:%M:%S %Y-%m-%d")
  echo -e "\n$time Filtering $sample ..."

  # filtration (-q 30: mapping quality>=30; -f 2: read mapped in proper pair)
  samtools view -b -q 30 -f 2 AlignmentRaw/${sample}/Aligned.sortedByCoord.out.bam > AlignmentFiltered/${sample}.bam

  # flagstat
  samtools flagstat AlignmentFiltered/${sample}.bam > QC/AlignmentFiltered/${sample}_flagstat.txt
done


#=====================================================================================================
# Step5: Count reads

for sample in ${Samples[@]};do
  time=$(date "+%H:%M:%S %Y-%m-%d")
  echo -e "\n$time Counting reads of $sample ..."

  featureCounts -p -a ${Gtf} -t exon -g gene_id -o ReadsCount/${sample}_readscount.txt -s 0 -T 12 AlignmentFiltered/${sample}.bam
  cut -f1,7 ReadsCount/${sample}_readscount.txt |sed "1i Gene\t${sample}" > ReadsCount/${sample}.tmp
done

# merge into reads count matrix
cd ReadsCount
paste Com1.tmp Com2.tmp Com3.tmp Com4.tmp D1.tmp D2.tmp D3.tmp D4.tmp MI1.tmp MI2.tmp MI3.tmp MI4.tmp |awk 'BEGIN {FS="\t";OFS="\t";} NR==1 || NR>=4 {print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24}' > readscount.txt
rm *tmp
cd ${WorkDir}


#=====================================================================================================
# Step6: Merge quality control results

multiqc QC/FastqRaw/*fastqc.zip QC/FastqTrimmed/*fastqc.zip AlignmentRaw/*/Log.final.out ReadsCount/*summary -o QC/multiqc

