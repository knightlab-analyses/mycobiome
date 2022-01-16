#!/bin/bash
#SBATCH --job-name=mlFungiTCGA8Cancers_WGS_RNA   # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=gdpoore@eng.ucsd.edu   # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --cpus-per-task=16            # Number of CPU cores per task
#SBATCH --mem=32gb                   # Job Memory
#SBATCH --time=24:00:00             # Time limit hrs:min:sec
#SBATCH --output=/home/gdpoore/cluster/logs/mlFungiTCGA8Cancers_WGS_RNA_%j.log    # Standard output and error log

work_dir=/home/gdpoore/projects/atlas_cancer_microbiome/AAA_Final_Analyses/AAA_DecontamV2/ML_mycobiome_tcga_8cancers_wgs_vs_rna
cd ${work_dir}

source activate tcgaAnalysisPythonR

Rscript S13B-ML-fungi-10k-rep1-tcga-8cancers-wgs-vs-rna-all-taxa-levels.R

