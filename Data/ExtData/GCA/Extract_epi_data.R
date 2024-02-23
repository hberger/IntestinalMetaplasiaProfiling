library(SeuratDisk)

data_folder = "."
Convert(source=file.path(data_folder, "epi_raw_counts02.h5ad"), dest = "h5seurat" )
d = LoadH5Seurat(file.path(data_folder,"epi_raw_counts02.h5seurat"))

d_healthy_intestine = subset(d, Diagnosis == "Healthy adult" & Region %in% c("LargeInt","SmallInt","REC","APD"))
SaveH5Seurat(d_healthy_intestine, filename=file.path(data_folder, "epi_raw_counts02_healthy_adult_intestine.h5seurat"), overwrite = TRUE)

