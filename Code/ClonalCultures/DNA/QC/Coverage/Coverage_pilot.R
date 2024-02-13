
setwd("/data_genome2/projects/MB044_gastric_clones/QC/DNA/Coverage/")
library(TEQC)
BAM_path="/data_genome2/projects/MB044_gastric_clones/mapping/DNA/BWA_MEM"
design_targets_file = "/data_genome1/References/AgilentSureSelect/Human_Exome_V6_UTR/S07604624_Covered_fixed_chroms.bed"

BAM_files = dir(path=BAM_path, pattern="*.bam$")

for (bf in BAM_files) {
  TEQCreport(sampleName=gsub("\\.bam","",bf),
  targetsName="Agilent SureSelect Human All Exon V6+UTR",
  referenceName="Human Genome GRCh37",
  destDir=paste("report",gsub("\\.bam","",bf),sep="/" ),
  reads=get.reads(file.path(BAM_path, bf)),
  targets=get.targets(design_targets_file, skip=0),
  genome="hg19")
}

report_folder_names = dir(path="report")

multiTEQCreport(singleReportDirs=paste("report",report_folder_names,sep="/"),
samplenames=report_folder_names,
projectName="MB044 Gastric clones",
targetsName="Agilent SureSelect Human All Exon V6+UTR",
referenceName="GRCh37",
destDir="multiTEQCreport")


