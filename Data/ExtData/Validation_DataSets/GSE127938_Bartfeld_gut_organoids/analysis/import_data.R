library(DESeq2)
library(data.table)
setwd("/data_genome2/public_data/GSE127938_Bartfeld_gut_organoids/analysis/")
mouse_data = fread("zcat ../Raw_counts/GSE127938_mouse_gastro_organoids_read_counts.txt.gz")
human_data = fread("zcat ../Raw_counts/GSE127938_human_gastro_organoids_read_counts.txt.gz")

human_mat = as.matrix(human_data[,3:ncol(human_data)])
rownames(human_mat) = human_data$V1
human_gene_anno = data.frame(GeneSymbol=human_data$GeneSymbol, row.names = as.character(human_data$V1), stringsAsFactors = F)

human_col_anno = data.frame(sampleID = colnames(human_data)[3:ncol(human_data)], stringsAsFactors = F)
tmp = strsplit(human_col_anno$sampleID, "_")
human_col_anno$species = "human"
tissue_trans = c("corp"="stomach corpus", "pyl"="stomach pylorus", "duod"="duodenum", "jej"="jejunum","ileum"="ileum","colon"="colon")
human_col_anno$tissue = tissue_trans[unlist(sapply(tmp,`[`,2))]
human_col_anno$replicate = unlist(sapply(tmp,`[`,4))
rownames(human_col_anno) = human_col_anno$sampleID


ds_human = DESeqDataSetFromMatrix( countData = human_mat, colData = human_col_anno, design = ~ tissue )

human_expr_mat = assay(vst(ds_human, blind=F))

save(human_expr_mat, human_col_anno, human_gene_anno, file="GSE127938_human_RNASeq.Rdata")

sel_genes = c("SIM2","MEIS1","MEIS2","CDX2","CDX1","VIL1","PDX1","IRX2","IRX5","GATA4", "GAST", "ATP4A", "PGC", "LEF1")
for (g in sel_genes) {
  v = human_expr_mat[rownames(subset(human_gene_anno, GeneSymbol==g)),]
  par(mfrow=c(1,2))
  barplot(v, las=2, main=g)
  boxplot(v ~ factor(human_col_anno[names(v),"tissue"], levels=c("stomach corpus","stomach pylorus","duodenum","jejunum","ileum","colon")) , las=2)
  par(mfrow=c(1,1))
}

mouse_mat = as.matrix(mouse_data[,3:ncol(mouse_data)])
rownames(mouse_mat) = mouse_data$V1
mouse_gene_anno = data.frame(GeneSymbol=mouse_data$GeneSymbol, row.names = as.character(mouse_data$V1), stringsAsFactors = F)

mouse_col_anno = data.frame(sampleID = colnames(mouse_data)[3:ncol(mouse_data)], stringsAsFactors = F)
tmp = strsplit(mouse_col_anno$sampleID, "_")
mouse_col_anno$species = "mouse"
tissue_trans = c("corp"="stomach corpus", "pyl"="stomach pylorus", "duod"="duodenum", "jej"="jejunum","ileum"="ileum","colon"="colon", "sto"="embryonal stomach", "si"="embryonal small intestine")
mouse_col_anno$tissue = tissue_trans[unlist(sapply(tmp,`[`,2))]
mouse_col_anno$replicate = unlist(sapply(tmp,`[`,4))
rownames(mouse_col_anno) = mouse_col_anno$sampleID


ds_mouse = DESeqDataSetFromMatrix( countData = mouse_mat, colData = mouse_col_anno, design = ~ tissue )

mouse_expr_mat = assay(vst(ds_mouse, blind=F))

save(mouse_expr_mat, mouse_col_anno, mouse_gene_anno, file="GSE127938_mouse_RNASeq.Rdata")

sel_genes = c("Sim2","Meis1","Meis2","Cdx2","Cdx1","Vil1","Pdx1","Irx2","Irx5","Gata4", "Gast", "Atp4a", "Pgc")
for (g in sel_genes) {
  v = mouse_expr_mat[rownames(subset(mouse_gene_anno, GeneSymbol==g)),]
  par(mfrow=c(1,2))
  barplot(v, las=2, main=g)
  boxplot(v ~ factor(mouse_col_anno[names(v),"tissue"], levels=c("stomach corpus","stomach pylorus","duodenum","jejunum","ileum","colon")) , las=2)
  par(mfrow=c(1,1))
}

###########################################
load("/data_genome2/public_data/GSE127938_Bartfeld_gut_organoids/analysis/GSE127938_human_RNASeq.Rdata")

sel_genes = c("PRSS1","PRSS2","CDX2","CDX1","MKI67", "PDX1", "HOXB3","HOXA2","HOXB5","ID3")

for (g in sel_genes) {
  v = human_expr_mat[rownames(subset(human_gene_anno, GeneSymbol==g)),]
  par(mfrow=c(1,2))
  barplot(v, las=2, main=g)
  boxplot(v ~ factor(human_col_anno[names(v),"tissue"], levels=c("stomach corpus","stomach pylorus","duodenum","jejunum","ileum","colon")) , las=2)
  par(mfrow=c(1,1))
}

