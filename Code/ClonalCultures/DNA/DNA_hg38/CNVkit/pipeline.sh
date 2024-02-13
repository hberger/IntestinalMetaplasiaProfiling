BAM_FOLDER=../../../mapping/DNA_hg38/BWA_MEM_hs_only_merged

CNVKIT_FOLDER=/data_genome1/SharedSoftware/CNVkit/download/cnvkit
REGIONS=/data_genome1/References/AgilentSureSelect/Human_Exome_V6_UTR/hg38/S07604624_Covered_fixed.bed

REF_FASTA=/data_genome1/References/Human/Sequences/Genome/GRCh38_no_alt_plus_hs38d1_Verily/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna

ACCESS_MAPPABLE_FILE=./access-10kb-mappable.hg38.bed

if [ ! -e $ACCESS_MAPPABLE_FILE ]; then 
  python ${CNVKIT_FOLDER}/cnvkit.py access $REF_FASTA -s 10000 -o $ACCESS_MAPPABLE_FILE
fi

all_chroms="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"


control="4447_M" # 591
samples=$(ls ${BAM_FOLDER}/4447_B_*.bam)
echo $samples

python ${CNVKIT_FOLDER}/cnvkit.py batch $samples --normal ${BAM_FOLDER}/${control}_*.bam \
    --targets $REGIONS --fasta $REF_FASTA \
    --access $ACCESS_MAPPABLE_FILE \
    --output-reference my_reference.cnn --output-dir results

control="4447_N" # 592
samples=$(ls ${BAM_FOLDER}/{4447_C,4447_D,4447_E,4447_F,4447_G,4447_H}_*.bam)
echo $samples

python ${CNVKIT_FOLDER}/cnvkit.py batch $samples --normal ${BAM_FOLDER}/${control}_*.bam \
    --targets $REGIONS --fasta $REF_FASTA \
    --access $ACCESS_MAPPABLE_FILE \
    --output-reference my_reference.cnn --output-dir results

control="4447_O" # 594
samples=$(ls ${BAM_FOLDER}/{4447_I,4447_J,4447_K,4447_L,4447_P,4447_Q}_*.bam)

python ${CNVKIT_FOLDER}/cnvkit.py batch $samples --normal ${BAM_FOLDER}/${control}_*.bam \
    --targets $REGIONS --fasta $REF_FASTA \
    --access $ACCESS_MAPPABLE_FILE \
    --output-reference my_reference.cnn --output-dir results

