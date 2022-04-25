#!/bin/bash
#SBATCH --job-name=mlFungiFullWgsRna   # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=gdpoore@eng.ucsd.edu   # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --cpus-per-task=32            # Number of CPU cores per task
#SBATCH --mem=64gb                   # Job Memory
#SBATCH --time=24:00:00             # Time limit hrs:min:sec
#SBATCH --output=/home/gdpoore/cluster/logs/mlFungiFullWgsRna_%j.log    # Standard output and error log

work_dir=/home/gdpoore/projects/atlas_cancer_microbiome/AA_Resubmission/ml_mycobiome_tcga_full_wgs_rna_vsnm_STonly
cd ${work_dir}

source activate tcgaAnalysisPythonR

Rscript S06R-ML-fungi-10k-rep1-tcga-full-wgs-rna-vsnm-STonly-VSNM-decontamV2.R

