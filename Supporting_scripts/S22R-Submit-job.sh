#!/bin/bash
#SBATCH --job-name=mlWISFungi   # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=gdpoore@eng.ucsd.edu   # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks-per-node=64            # Number of CPU cores per task
#SBATCH --mem=64gb                   # Job Memory
#SBATCH --time=100:00:00             # Time limit hrs:min:sec
#SBATCH --output=/home/gdpoore/cluster/logs/mlWISFungi_%j.log    # Standard output and error log

work_dir=/home/gdpoore/projects/atlas_cancer_microbiome/AA_Resubmission/ml_tcga_wis_intersected
cd ${work_dir}

source activate r4

Rscript S22R-ML-fungi-wis-vs-nonwis.R

