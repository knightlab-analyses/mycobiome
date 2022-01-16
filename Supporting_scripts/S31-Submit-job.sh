#!/bin/bash
#SBATCH --job-name=mlBDNEarlyVSNM   # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=gdpoore@eng.ucsd.edu   # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --cpus-per-task=32            # Number of CPU cores per task
#SBATCH --mem=64gb                   # Job Memory
#SBATCH --time=24:00:00             # Time limit hrs:min:sec
#SBATCH --output=/home/gdpoore/cluster/logs/mlBDNEarlyVSNM_%j.log    # Standard output and error log

work_dir=/home/gdpoore/projects/atlas_cancer_microbiome/AAA_Final_Analyses/AAA_DecontamV2/ML_BDN_early_raw
cd ${work_dir}

source activate tcgaAnalysisPythonR

Rscript S31-ML-fungi-10k-rep1-tcga-bdn-early-stage-by-seqcenter.R

