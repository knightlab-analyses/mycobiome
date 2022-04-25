#!/bin/bash
#SBATCH --job-name=mlEukDetectHMP   # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=gdpoore@eng.ucsd.edu   # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --cpus-per-task=16            # Number of CPU cores per task
#SBATCH --mem=32gb                   # Job Memory
#SBATCH --time=48:00:00             # Time limit hrs:min:sec
#SBATCH --output=/home/gdpoore/cluster/logs/mlEukDetectHMP_%j.log    # Standard output and error log

work_dir=/home/gdpoore/projects/atlas_cancer_microbiome/AA_Resubmission/ml_tcga_seq_center_eukdetect_hmp_taxa_levels
cd ${work_dir}

source activate tcgaAnalysisPythonR

Rscript S17R-ML-fungi-10k-rep1-tcga-eukdetect-hmp-by-seqcenter-taxa-levels.R

