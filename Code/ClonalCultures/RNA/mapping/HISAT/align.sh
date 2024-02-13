#!/bin/bash
#ulimit -v 30258917744
HISAT_FOLDER=/data_genome1/SharedSoftware/HISAT/download/hisat2-2.1.0
GENOME_DIR=/data_genome1/References/GenomeIndices/Human/HISAT2/hg38/genome 
SPLICE_SITES=/data_genome1/References/Human/Annotations/hg38/Ensembl/splice_sites.txt

export PATH=$PATH:$HISAT_FOLDER

function join { local IFS="$1"; shift; echo "$*"; }

TMP_FOLDER=./tmp

mkdir $TMP_FOLDER

cat ../../../../metadata/phenotyping/v2_2020_01/all_samples_aggregated.txt | while read sample_line
do # 1st pass
  [ -z "$sample_line" ] && continue
  set $sample_line
  sample=$(echo $1)
  reads1=$(echo $2)
#  reads2=$(echo $3)
#  $STARexec --genomeDir $GENOME_DIR --readFilesIn $reads1 --readFilesCommand zcat --outFileNamePrefix ${sample}_ --runThreadN 8 --outSAMstrandField intronMotif --genomeLoad NoSharedMemory --twopassMode Basic

  $HISAT_FOLDER/hisat2 -p 12 --known-splicesite-infile $SPLICE_SITES  -x $GENOME_DIR -U $reads1 -S ${sample}_Aligned.out.sam

  /data_genome1/SharedSoftware/samtools_htslib/samtools view -u -b ${sample}_Aligned.out.sam | /data_genome1/SharedSoftware/samtools_htslib/samtools sort -@ 4 -T ${TMP_FOLDER}/CHUNK -O bam -o ${sample}_sorted.bam
  /data_genome1/SharedSoftware/samtools_htslib/samtools index -b ${sample}_sorted.bam

  rm ${sample}_Aligned.out.sam

done

