setwd("~/Work/UKSH/IM/IntestinalMetaplasiaProfiling/Code/ClonalCultures/RNA/DGE_analysis_final/")

library(org.Hs.eg.db)
library(AnnotationDbi)
library(data.table)
library(readxl)
  
count_data_hisat = fread("../../../../Data/Processed/ClonalCultures/RNA/featureCount/Ensembl.98/ESGIv2_counts.txt.gz")
sample_cols = grepl(".bam$", colnames(count_data_hisat), perl=T)

count_mat = as.matrix(count_data_hisat[, sample_cols, with=F])
rownames(count_mat) = count_data_hisat$Geneid
snames = gsub("_sorted\\.bam", "", unlist(sapply(strsplit(colnames(count_mat) |> basename(), "/"), function(x) x[length(x)] )) )
colnames(count_mat) = snames

gene_anno = as.data.frame(count_data_hisat[, !sample_cols, with=F])


# this was run in 2020
# library(biomaRt)
# ensembl90 <- useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', version = 98)
#  
# gg = getBM(attributes = c("entrezgene","hgnc_symbol","ensembl_gene_id", "description"), filter = c("chromosome_name", "ensembl_gene_id"), values=list(chromosome_name = c(as.character(1:22), "X","Y"), ensembl_gene_id = rownames(count_mat)), mart=ensembl90 ) 
#  
# fwrite(gg, file="../../../../Data/ExtData/Ensembl/Ensembl_v98_human_gene_annotation.txt", sep="\t", quote=F)

ens_anno_orig = fread("../../../../Data/ExtData/Ensembl/Ensembl_v98_human_gene_annotation.txt")
setnames(ens_anno_orig, colnames(ens_anno_orig), c("EntrezID","GeneSymbol", "Geneid", "GeneName"))
# as usuaö there are 1:n relationships to GeneSymbol, EntrezID, GeneName - fix those to 1:1
collapse_fun <- function(x) paste(unique(x[!is.na(x)]), collapse="," )
ens_anno_agg = ens_anno_orig[, .(EntrezID=collapse_fun(EntrezID), GeneSymbol=collapse_fun(GeneSymbol), GeneName=collapse_fun(GeneName)), by = "Geneid" ]

gene_anno = merge(gene_anno, ens_anno_agg, by="Geneid", all.x=T, sort=F)
rownames(gene_anno) = gene_anno$Geneid

# FIXME consider reading total counts (including unmapped/unassigned) for normalization
fpkm_mat = sweep(sweep(count_mat, 2, apply(count_mat, 2, sum)/1e6, "/"), 1, gene_anno[rownames(count_mat), ]$Length/1e3, "/")

exp_design = as.data.frame(read_excel("../../../../Data/Metadata/ClonalCultures/old/phenotyping/SampleDescription.xlsx", sheet = 1), stringsAsFactors=FALSE)
colnames(exp_design) = make.names(colnames(exp_design))
exp_design = subset(exp_design, !is.na(Library.number))
exp_design$sampleID = as.character(exp_design$Library.number)
rownames(exp_design) = exp_design$sampleID
ed = exp_design


save(count_mat, fpkm_mat, gene_anno, ed, file="../../../../Data/Processed/ClonalCultures/RNA/HISAT_featureCount_hg38_gene_expression.Rdata", compress=T)

