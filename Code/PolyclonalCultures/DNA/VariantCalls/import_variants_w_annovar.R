library(VariantAnnotation)
library(data.table)

setwd("Code/PolyclonalCultures/DNA/VariantCalls/")

all_samples_folder = "../../../../Data/Processed/PolyclonalCultures/WES/VariantCalls/"

all_vcf_raw = list()
# annotated variants, all samples
all_vcf_raw[["ZZ_IM"]] = readVcf(file.path(all_samples_folder, "IManno.vcf"))

all_samples = list()
# all_samples[["ZZ_IM"]] = c("IM10A-W-R","IM10",
#                            "IM11A-W-R","IM11",
#                            "IM13A-W-R","IM13A",
#                            "IM19A-W-R","non19C-W-R",
#                            "IM20A-W-R","non20C-W-R",
#                            "IM30A-W-R","IM30",
#                            "IM43AC-W-R","IM43-blood",
#                            "IM51A-W-R","IM51",
#                            "IM52-W-R","IM52",
#                            "IM54A-W-R","non54C-W-R")


all_samples[["ZZ_IM"]] = c("IM10A-W-R",
                           "IM11A-W-R",
                           "IM13A-W-R",
                           #"IM19A-W-R", # target or control DNA mixup; most variants (n = 16109) are at known SNP sites
                           "IM20A-W-R",
                           "IM30A-W-R",
                           "IM43AC-W-R",
                           "IM51A-W-R",
                           "IM52-W-R",
                           "IM54A-W-R")


# Annotations from SnpEff are per matching transcript, i.e. there can be more than one annotations per line. 
# This must be processed appropriately, e.g. by the script vcf_to_ts_csv_SnpEff_plus_filters.py
#ann_header_info = as.data.frame(info(header(vcf_raw)))["ANN","Description"]
#aa = regexpr("'[^\\']*'",ann_header_info, perl=T)
#ann_cols = make.names(unlist(strsplit(substr(ann_header_info, aa+1, aa + attr(aa, "match.length")-2), " \\| ")))
#tmp = info(vcf_raw)[["ANN"]]

# import SnpEff annotated variants
get_alt_allele_depth <- function(var_allele, a, c, g, t) {
  switch(var_allele, 
         A=a,
         C=c,
         G=g,
         T=t,
         NA
  )
}

DNAstringlist2char = function(x) sapply(lapply(x, as.character), paste, collapse = ",") 
DNAstringlist2maxlen = function(x) sapply(lapply(x, as.character), function(y) max(nchar(y))) 

first = 1
all_variants = c()
for (expid in names(all_vcf_raw)) {

  vcf_raw = all_vcf_raw[[expid]]

  r_tmp = rowRanges(vcf_raw)
  
  r_tmp$multiallelic = unlist(lapply(r_tmp$ALT, length))>1
  r_tmp$alt_allele = unlist(lapply(r_tmp$ALT, paste, collapse=","))

  maxlen_alt = DNAstringlist2maxlen(rowRanges(vcf_raw)$ALT)
  maxlen_ref = DNAstringlist2maxlen(rowRanges(vcf_raw)$REF)
  
  vartype = ifelse(maxlen_ref >1 | maxlen_alt > 1, "Indel", "SNV")
  r_tmp$vartype = vartype
  
  dp = geno(vcf_raw)[["DP"]]
  ad = geno(vcf_raw)[["AD"]]

  for (s in all_samples[[expid]]) {
    r = r_tmp
    r$sampleID = s
    r$DP = dp[,s]

    r$AC = unlist(lapply((ad[,s]),`[`,2))

    # remove variants not detected in this sample
    r = r[!is.na(r$AC)]
    
    if (first) {
      all_variants = r
      first=0
    } else {
      all_variants = c(all_variants, r)
    }
  }
}

# Variant Annotation

# import Annovar annotations

annovar_output_folder = all_samples_folder 
annovar_output_file = "calls_merged_ANNOVAR.hg38_multianno.txt"
annovar_data = read.table(file.path(annovar_output_folder,annovar_output_file), sep="\t", fill=T, header=F, skip = 1, stringsAsFactors = F, na.strings = c("NA", "."), quote = "")
hh = read.table(file.path(annovar_output_folder,annovar_output_file), sep="\t", fill=T, header=F, skip = 0, nrows = 1, stringsAsFactors = F)
colnames(annovar_data)[1:ncol(hh)] = hh[1,]
rm(hh)
annovar_data$ID = with(annovar_data, paste0(Chr,"_",Start,"_",Ref,"_",Alt))

# unfortunately annovar normalizes variants, channging start/ref/alt for indels. 
# We have to normalized to un-normalized variant descriptions in order to match to SnpEff
annovar_ids = fread(file.path(annovar_output_folder, "calls_merged_ANNOVAR.avinput"))
annovar_ids$ID_normalized = with(annovar_ids, paste0(V1,"_",V2,"_",V4,"_",V5) )
annovar_ids$ID_unnormalized = with(annovar_ids, paste0(V9,"_",V10,"_",V12,"_",V13) )
annovar_data = merge(annovar_data, annovar_ids[, c("ID_normalized","ID_unnormalized")], by.x="ID", by.y="ID_normalized", all.x=T, sort=F)

# We now take SnpEff annotations already preprocessed by a custom python script
all_anno = list()
tmp = fread(file.path(all_samples_folder,"Variant_annotation_ZZ_IM_hg38_Mutect2.txt"), sep="\t", header=T, stringsAsFactors = F)
tmp$ID = with(tmp, paste0(CHROM,"_",POS,"_",REF,"_",ALT))
tmp2 = merge(tmp, annovar_data, by.x="ID", by.y="ID_unnormalized", all.x=T, sort=F)
tmp2[["ID.y"]] <- NULL

all_anno[["ZZ_IM"]] = tmp2

impact_order = c("HIGH"=4,"MODERATE"=3,"LOW"=2,"MODIFIER"=1)

all_anno_df = unique(do.call(rbind, all_anno))
all_anno_df[,IMPACT_ORDER := impact_order[Impact]]
all_anno_df[,MAX_IMPACT := max(IMPACT_ORDER), by=c("CHROM","POS","REF","ALT")]

anno_unique_sites = all_anno_df[,.(Gene=list(Gene), Impact=list(Impact), AnnotationType=list(AnnotationType), dbSNP=unique(dbSNP) ),by=c("CHROM","POS","REF","ALT")]

anno_unique_sites_highest_impact = all_anno_df[IMPACT_ORDER==MAX_IMPACT,][,.(Gene=list(Gene), ENSG=list(ENSG), HighestImpact=unique(Impact), AnnotationType=list(AnnotationType), dbSNP=unique(dbSNP) ),by=c("CHROM","POS","REF","ALT")]
anno_unique_sites_highest_impact$ID = with(anno_unique_sites_highest_impact, paste0(CHROM,"_",POS,"_",REF,"_",ALT))
anno_unique_sites_highest_impact = merge(anno_unique_sites_highest_impact, annovar_data, by.x="ID", by.y="ID_unnormalized", all.x=T, sort=F)
anno_unique_sites_highest_impact$SnpEff_annovar_gene_ovlp =  mapply(function(x,y) {length(intersect(x,y))} , anno_unique_sites_highest_impact$Gene, anno_unique_sites_highest_impact$Gene.refGene)

#zz = subset(anno_unique_sites_highest_impact_w_annovar, Gene != Gene.refGene)
zz = subset(anno_unique_sites_highest_impact, SnpEff_annovar_gene_ovlp==0)

# Discrepant variant annotations (distinct gene symbol)
table(zz$SnpEff_annovar_gene_ovlp==0, zz$HighestImpact)
# many of the discrepancies can be explained by different gene models (ensembl vs. RefGene) or different symbols for the same gene

save(all_variants, all_samples, anno_unique_sites, anno_unique_sites_highest_impact, file=file.path(all_samples_folder, "ZZ_IM_Exome_Mutect2_hg38_Calls.Rdata"))

