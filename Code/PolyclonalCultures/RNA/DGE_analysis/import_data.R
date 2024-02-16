setwd("~/Work/UKSH/IM/IntestinalMetaplasiaProfiling/Code/PolyclonalCultures/RNA/DGE_analysis")

library(org.Hs.eg.db)
library(AnnotationDbi)
library(data.table)
library(readxl)
  
count_data_hisat = fread("../../../../Data/Processed/PolyclonalCultures/RNA/featureCount_Ens98/ESGIv2_counts_PC_ZZ.txt.gz")
sample_cols = grepl(".bam$", colnames(count_data_hisat), perl=T)

count_mat = as.matrix(count_data_hisat[, sample_cols, with=F])
rownames(count_mat) = count_data_hisat$Geneid
snames = gsub("_sorted\\.bam", "", unlist(sapply(strsplit(colnames(count_mat) |> basename(), "/"), function(x) x[length(x)] )) )
colnames(count_mat) = snames

gene_anno = as.data.frame(count_data_hisat[, !sample_cols, with=F])


# this was run in 2020
# library(biomaRt)
# ensembl98 <- useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', version = 98)
#  
# gg = getBM(attributes = c("entrezgene","hgnc_symbol","ensembl_gene_id", "description"), filter = c("chromosome_name", "ensembl_gene_id"), values=list(chromosome_name = c(as.character(1:22), "X","Y"), ensembl_gene_id = rownames(count_mat)), mart=ensembl98 ) 
#  
# fwrite(gg, file="../../../../Data/ExtData/Ensembl/Ensembl_v98_human_gene_annotation.txt", sep="\t", quote=F)

ens_anno_orig = fread("../../../../Data/ExtData/Ensembl/Ensembl_v98_human_gene_annotation.txt")
setnames(ens_anno_orig, colnames(ens_anno_orig), c("EntrezID","GeneSymbol", "Geneid", "GeneName"))
# as usual there are 1:n relationships to GeneSymbol, EntrezID, GeneName - fix those to 1:1
collapse_fun <- function(x) paste(unique(x[!is.na(x)]), collapse="," )
ens_anno_agg = ens_anno_orig[, .(EntrezID=collapse_fun(EntrezID), GeneSymbol=collapse_fun(GeneSymbol), GeneName=collapse_fun(GeneName)), by = "Geneid" ]

gene_anno = merge(gene_anno, ens_anno_agg, by="Geneid", all.x=T, sort=F)
rownames(gene_anno) = gene_anno$Geneid

# FIXME consider reading total counts (including unmapped/unassigned) for normalization
fpkm_mat = sweep(sweep(count_mat, 2, apply(count_mat, 2, sum)/1e6, "/"), 1, gene_anno[rownames(count_mat), ]$Length/1e3, "/")

exp_design = as.data.frame(read_xlsx("../../../../Data/Metadata/PolyclonalCultures/SampleDescription_phenotyping.xlsx", sheet = 1), stringsAsFactors=FALSE)
exp_design$Group = ifelse(grepl("^NON IM",exp_design$Remark), "Non-IM", ifelse(grepl("^IM", exp_design$Remark),"IM", "Control"))
exp_design$Patient = gsub("-|\\s","", gsub("^IM|^NON IM","", gsub("A$|AC$|C$","", exp_design$Remark)))
rownames(exp_design) = exp_design$sampleID
ed = exp_design

save(count_mat, fpkm_mat, gene_anno, ed, file="../../../../Data/Processed/PolyclonalCultures/RNA/HISAT_featureCount_hg38_gene_expression.Rdata", compress=T)

