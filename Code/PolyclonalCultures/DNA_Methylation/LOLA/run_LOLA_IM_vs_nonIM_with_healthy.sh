#!/bin/bash
#
#SBATCH --job-name=CellChat_complete_NAFLD
#SBATCH --output=slurm_output_complete_NAFLD
#
#SBATCH --ntasks=12
#SBATCH --nodes=1
#SBATCH --time=80:00:00
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hilmar.berger@charite.de

r_file=run_LOLA_IM_vs_nonIM_with_healthy.R
R=~/shared/R/4.2.2/bin/R

(cat <<"EOF"
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2")
library(knitr)
rmarkdown::render("LOLA_IM_vs_nonIM_with_healthy.Rmd")
quit("no")
EOF
) > $r_file


$R CMD BATCH $r_file

