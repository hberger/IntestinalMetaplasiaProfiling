ESGIv2_calls=../Mutect2/called_variants/SnpEff/ESGI_v2_GRCh38_WXS_calls_merged_SnpEff.vcf.gz

python ./vcf_to_ts_anno_SnpEff_plus_filters.py -o Variant_annotation_ESGIv2_hg38_Mutect2.txt $ESGIv2_calls

