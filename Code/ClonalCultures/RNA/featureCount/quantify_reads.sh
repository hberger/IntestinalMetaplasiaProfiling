#!/bin/bash
ulimit -v 30258917744
#INDEX_FOLDER=/data_genome1/References/GenomeIndices/Human/Transcriptome_index/Gencode/v25lift37_Salmon_v1.1.0/GencodeHv25lift37
GENE_MODEL=/data_genome1/References/Human/Annotations/hg38/Ensembl/Homo_sapiens.GRCh38.98.chr.gtf

function join { local IFS="$1"; shift; echo "$*"; }

FEATURECOUNT_EXEC=/data_genome1/SharedSoftware/Subread/bin/featureCounts
TMP_FOLDER=./tmp
#GENE_MODEL=/data_genome1/References/Human/Annotations/hg19/Gencode/v25/gencode.v25lift37.annotation.gtf
BAM_FOLDER=../../../mapping/RNA/phenotyping/HISAT

mkdir $TMP_FOLDER

#cat ../../../metadata/phenotyping/v2_2020_01/all_samples_aggregated.txt | while read sample_line
#do # 1st pass
#  [ -z "$sample_line" ] && continue
#  set $sample_line
#  sample=$(echo $1)
#  bam_files=$(ls ${BAM_FOLDER}/${sample}_*.bam )
#  reads1=$(echo $2 | sed -e 's/,/ /g')
#  reads2=$(echo $3)

  #$SALMON_EXEC quant -i $INDEX_FOLDER -l A -r $reads1 -o $sample -g $GENE_MODEL --validateMappings
  $FEATURECOUNT_EXEC -M --primary -s 2 -T 8 -t exon -g gene_id -a $GENE_MODEL -o ESGIv2_counts.txt $BAM_FOLDER/*_sorted.bam

#done

