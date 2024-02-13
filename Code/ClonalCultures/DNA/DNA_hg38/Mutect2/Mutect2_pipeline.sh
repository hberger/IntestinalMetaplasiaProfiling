BAM_FOLDER=../GATK

GATK=/data_genome1/SharedSoftware/GATK/gatk-4.1.9.0/gatk

#REF_FASTA=/data_genome1/References/MusMusculus/Sequences/Genome/GRCm39/GCA_000001635.9_GRCm39_full_analysis_set.fna
REF_FASTA=/data_genome1/References/Human/Sequences/Genome/GRCh38_no_alt_plus_hs38d1_Verily/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna

REGIONS=/data_genome1/References/AgilentSureSelect/Human_Exome_V6_UTR/hg38/S07604624_Covered_fixed.bed
GNOMAD=/data_genome1/SharedSoftware/GATK/resources_Mutect2/af-only-gnomad.hg38.vcf.gz

function run_mutect2 { sample=$1;
           # --germline-resource germline_m2.vcf.gz does not work without AF field :(
	   $GATK Mutect2 \
           -R $REF_FASTA \
           -I $(ls ${BAM_FOLDER}/${sample}_merged_dedup.bam) \
           -tumor $sample \
           -I $(ls ${BAM_FOLDER}/${control}_merged_dedup.bam) \
           -normal $control \
	   --disable-adaptive-pruning \
	   --f1r2-tar-gz f1r2.tar.gz \
           --germline-resource $GNOMAD \
           -L $REGIONS \
           -O ${sample}.vcf.gz &> ${sample}.log; 


	   $GATK LearnReadOrientationModel -I f1r2.tar.gz -O read-orientation-model.tar.gz

	   #$GATK GetPileupSummaries \
    	   #-I tumor.bam \
           #-V chr17_small_exac_common_3_grch38.vcf.gz \
           #-L chr17_small_exac_common_3_grch38.vcf.gz \
           #-O getpileupsummaries.table

           #$GATK CalculateContamination \
           #-I getpileupsummaries.table \
           #-tumor-segmentation segments.table \
           #-O calculatecontamination.table

           $GATK FilterMutectCalls -V ${sample}.vcf.gz -R $REF_FASTA --ob-priors read-orientation-model.tar.gz -O ${sample}_filtered.vcf.gz

}

export -f run_mutect2
export DATA_FOLDER
export REF_FASTA
export BAM_FOLDER
export GATK
export GNOMAD
export REGIONS


control="4447_M" # 591
samples="4447_B"

for sample in $samples; do
    #export sample
    export control 

    #echo $sample | sed -e 's/ /\n/g' | parallel -t -j 12 run_mutect2
    run_mutect2 $sample
    bcftools view -f .,PASS -o ${sample}_.pass.vcf.gz -O z ${sample}_filtered.vcf.gz
    tabix -p vcf ${sample}_pass.vcf.gz
done

control="4447_N" # 592
samples="4447_C 4447_D 4447_E 4447_F 4447_G 4447_H"

for sample in $samples; do
    #export sample
    export control 

    #echo $sample | sed -e 's/ /\n/g' | parallel -t -j 12 run_mutect2
    run_mutect2 $sample
    bcftools view -f .,PASS -o ${sample}_pass.vcf.gz -O z ${sample}_filtered.vcf.gz
    tabix -p vcf ${sample}_pass.vcf.gz
done



control="4447_O" # 594
samples="4447_I 4447_J 4447_K 4447_L 4447_P 4447_Q"


for sample in $samples; do
    #export sample
    export control 

    #echo $sample | sed -e 's/ /\n/g' | parallel -t -j 12 run_mutect2
    run_mutect2 $sample
    bcftools view -f .,PASS -o ${sample}_pass.vcf.gz -O z ${sample}_filtered.vcf.gz
    tabix -p vcf ${sample}_pass.vcf.gz
done


