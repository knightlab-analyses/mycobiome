#!/bin/bash
#SBATCH --job-name=fungi_cov   # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=gdpoore@eng.ucsd.edu   # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --cpus-per-task=2            # Number of CPU cores per task
#SBATCH --mem=64gb                   # Job Memory
#SBATCH --time=24:00:00             # Time limit hrs:min:sec
#SBATCH --output=/home/gdpoore/cluster/logs/fungi_cov_%j.log    # Standard output and error log

in_dir=/projects/tcga-data/mycobiome_coverage_analysis/extracted_sam_files_fungi_hiseq_only
zebra_dir=/projects/tcga-data/mycobiome_coverage_analysis/zebra_filter
out_dir=/projects/tcga-data/mycobiome_coverage_analysis
cd ${working_dir}

source activate qiime2-2020.8

python ${zebra_dir}/calculate_coverages.py \
 -i ${in_dir}/ \
 -o ${out_dir}/fungi_filt_updated_29Sep21_coverage_all_wgs_and_rna_output_hiseq_only.tsv \
 -d ${zebra_dir}/databases/rep200/rep200_metadata.tsv