#!/bin/bash
#SBATCH --job-name=mlWIS   # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=gdpoore@eng.ucsd.edu   # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --cpus-per-task=32            # Number of CPU cores per task
#SBATCH --mem=64gb                   # Job Memory
#SBATCH --time=8:00:00             # Time limit hrs:min:sec
#SBATCH --output=/home/gdpoore/cluster/logs/mlWIS_%j.log    # Standard output and error log

work_dir=/home/gdpoore/projects/atlas_cancer_microbiome/AA_Resubmission/ml_wis_samples
cd ${work_dir}

source activate tcgaAnalysisPythonR

Rscript S26R-ML-WIS-fungi-bacteria.R

