#!/bin/bash
#SBATCH --job-name=PVCA   # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=gdpoore@eng.ucsd.edu   # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --cpus-per-task=16            # Number of CPU cores per task
#SBATCH --mem=64gb                   # Job Memory
#SBATCH --time=24:00:00             # Time limit hrs:min:sec
#SBATCH --output=/home/gdpoore/cluster/logs/PVCA_%j.log    # Standard output and error log

work_dir=/home/gdpoore/projects/atlas_cancer_microbiome/AAA_Final_Analyses/PVCA_mycobiome_taxa_levels
cd ${work_dir}

source activate tcgaAnalysisPythonR

Rscript S08-Run-pvca-all-taxa-levels.R

