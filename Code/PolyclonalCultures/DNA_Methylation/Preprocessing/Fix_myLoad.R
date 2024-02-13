# Run this from the main folder of the repository
load("./Data/Processed/PolyclonalCultures/DNAm/myLoad_IM_nonIM_and_healthy_stomach_combined.RData")
# for some reason the methylation object has less rows than expected for EPIC arrays while the rgSet object is correct - maybe the wrong annotation was available
# let's repeat some of the steps done my champ.load to get the full methylation data from EPIC arrays

detPcut <- 0.01
SampleCutoff <- 0.1
ProbeCutoff <- 0
beadCutoff <- 0.05

mset <- preprocessRaw(myLoad$rgSet)
detP <- detectionP(myLoad$rgSet)
tmp = getBeta(mset, "Illumina")
tmp[detP >= detPcut] <- NA
message("The fraction of failed positions per sample\n \n            (You may need to delete samples with high proportion of failed probes\n): ")
numfail <- matrix(colMeans(is.na(tmp)))
rownames(numfail) <- colnames(detP)
colnames(numfail) <- "Failed CpG Fraction."
print(numfail)

RemainSample <- which(numfail < SampleCutoff)
if (any(numfail >= SampleCutoff)) 
  message("The detSamplecut parameter is : ", SampleCutoff, 
          "\nSamples : ", paste(rownames(numfail)[which(numfail >= 
                                                          SampleCutoff)], collapse = ","), " will be deleted.\n", 
          "There are ", length(RemainSample), " samples left for analysis.\n")
rgSet <- myLoad$rgSet[, RemainSample]
detP <- detP[, RemainSample]
mset <- mset[, RemainSample]
pd <- myLoad$pd[RemainSample, ]
tmp <- tmp[, RemainSample]

# Filtering
mset.f = mset[rowSums(is.na(tmp)) <= ProbeCutoff * 
                ncol(detP), ]
message("Filtering probes with a detection p-value above ", 
        detPcut, " in one or more samples has removed ", 
        dim(mset)[1] - dim(mset.f)[1], " probes from the analysis. If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples.")
mset = mset.f
tmp <- tmp[rowSums(is.na(tmp)) <= ProbeCutoff * ncol(detP), ]


mybeadcount <- function(x) {
  nb <- getNBeads(x)
  typeIadd <- getProbeInfo(x, type = "I")
  typeImatchA <- match(typeIadd$AddressA, rownames(nb))
  typeImatchB <- match(typeIadd$AddressB, rownames(nb))
  typeIIadd <- getProbeInfo(x, type = "II")
  typeIImatch <- match(typeIIadd$Address, rownames(nb))
  nbcg <- nb
  locusNames <- getManifestInfo(x, "locusNames")
  bc_temp <- matrix(NA_real_, ncol = ncol(x), nrow = length(locusNames), 
                    dimnames = list(locusNames, sampleNames(x)))
  TypeII.Name <- getProbeInfo(x, type = "II")$Name
  bc_temp[TypeII.Name, ] <- nbcg[getProbeInfo(x, type = "II")$AddressA, 
  ]
  TypeI <- getProbeInfo(x, type = "I")
  bcB <- bc_temp
  bcA <- bc_temp
  bcB[TypeI$Name, ] <- nbcg[TypeI$AddressB, ]
  bcA[TypeI$Name, ] <- nbcg[TypeI$AddressA, ]
  bcB3 <- which(bcB < 3)
  bcA3 <- which(bcA < 3)
  bcA2 <- bcA
  bcB2 <- bcB
  bcA2[bcA3] <- NA
  bcA2[bcB3] <- NA
  bc <- data.frame(bcA2)
  bc
}

# Filter bead counts
bc = mybeadcount(rgSet)
bc2 = bc[rowSums(is.na(bc)) < beadCutoff * (ncol(bc)), ]
mset.f2 = mset[featureNames(mset) %in% row.names(bc2), ]
tmp <- tmp[rownames(tmp) %in% row.names(bc2), ]
message("Filtering probes with a beadcount <3 in at least ", 
        beadCutoff * 100, "% of samples, has removed ", 
        dim(mset)[1] - dim(mset.f2)[1], " from the analysis.")
mset = mset.f2

# Filter non-CG
mset.f2 = dropMethylationLoci(mset, dropCH = T)
tmp <- tmp[rownames(tmp) %in% featureNames(mset.f2), ]
message("Filtering non-cg probes, has removed ", 
        dim(mset)[1] - dim(mset.f2)[1], " from the analysis.")
mset <- mset.f2

# SNP filtering
message("Using general EPIC SNP list for filtering.")
data(EPIC.manifest.hg19)
maskname <- rownames(EPIC.manifest.hg19)[which(EPIC.manifest.hg19$MASK_general == TRUE)]
mset.f2 = mset[!featureNames(mset) %in% maskname, ]
tmp <- tmp[!rownames(tmp) %in% maskname, ]
message("Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper, 2016, has removed ", 
        dim(mset)[1] - dim(mset.f2)[1], " from the analysis.")
mset = mset.f2

# Filter Multi-hit
data(multi.hit)
mset.f2 = mset[!featureNames(mset) %in% multi.hit$TargetID, ]
tmp <- tmp[!rownames(tmp) %in% multi.hit$TargetID, ]
message("Filtering probes that align to multiple locations as identified in Nordlund et al, has removed ", 
        dim(mset)[1] - dim(mset.f2)[1], " from the analysis.")
mset = mset.f2

# Filter X/Y probes
data(probe.features.epic)
autosomes = probe.features[!probe.features$CHR %in% 
                             c("X", "Y"), ]
mset.f2 = mset[featureNames(mset) %in% row.names(autosomes), ]
tmp <- tmp[rownames(tmp) %in% row.names(autosomes), ]
message("Filtering probes on the X or Y chromosome has removed ", 
        dim(mset)[1] - dim(mset.f2)[1], " from the analysis.")
mset = mset.f2

beta.raw <- tmp
intensity <- minfi::getMeth(mset) + minfi::getUnmeth(mset)
detP <- detP[which(row.names(detP) %in% row.names(beta.raw)), ]

if (min(beta.raw, na.rm = TRUE) <= 0) 
  beta.raw[beta.raw <= 0] <- min(beta.raw[beta.raw > 
                                            0])
message("Zeros in your dataset have been replaced with smallest positive value.\n")
if (max(beta.raw, na.rm = TRUE) >= 0) 
  beta.raw[beta.raw >= 1] <- max(beta.raw[beta.raw < 
                                            1])
message("One in your dataset have been replaced with largest value below 1.\n")

message("The analysis will be proceed with ", dim(beta.raw)[1], 
        " probes and ", dim(beta.raw)[2], " samples.\n")
message("Current Data Set contains ", sum(is.na(beta.raw)), 
        " NA in ", "[Beta]", " Matrix.\n")

myLoad <- list(mset = mset, rgSet = rgSet, pd = pd, intensity = intensity, 
               beta = beta.raw, detP = detP)

save(myLoad, file="./Data/Processed/PolyclonalCultures/DNAm/myLoad_IM_nonIM_and_healthy_stomach_combined_fixed.RData")