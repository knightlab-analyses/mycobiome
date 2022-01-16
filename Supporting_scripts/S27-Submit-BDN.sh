#!/bin/bash
#SBATCH --job-name=mlBDNFungiVsBact   # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=gdpoore@eng.ucsd.edu   # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --cpus-per-task=16            # Number of CPU cores per task
#SBATCH --mem=32gb                   # Job Memory
#SBATCH --time=72:00:00             # Time limit hrs:min:sec
#SBATCH --output=/home/gdpoore/cluster/logs/mlBDNFungiVsBact_%j.log    # Standard output and error log

work_dir=/home/gdpoore/projects/atlas_cancer_microbiome/AAA_Final_Analyses/AAA_DecontamV2/ML_fungi_vs_bacteria_v2
cd ${work_dir}

source activate tcgaAnalysisPythonR

Rscript S27-ML-fungi-vs-bacteria-number-features-per-cancer-type-ranked-BDN.R

