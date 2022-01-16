#!/bin/bash
#SBATCH --job-name=mlWISFungiVsBact   # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=gdpoore@eng.ucsd.edu   # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --cpus-per-task=64            # Number of CPU cores per task
#SBATCH --mem=64gb                   # Job Memory
#SBATCH --time=48:00:00             # Time limit hrs:min:sec
#SBATCH --output=/home/gdpoore/cluster/logs/mlWISFungiVsBact_%j.log    # Standard output and error log

work_dir=/home/gdpoore/projects/atlas_cancer_microbiome/AAA_Final_Analyses/AAA_DecontamV2/ML_WIS_fungi_vs_bacteria
cd ${work_dir}

source activate r4

Rscript S24-ML-WIS-fungi-vs-bacteria-number-features.R

