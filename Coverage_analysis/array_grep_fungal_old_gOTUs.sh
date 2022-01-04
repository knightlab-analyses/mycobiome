#!/bin/bash
#SBATCH --job-name=fungi_extract   # Job name
#SBATCH --mail-type=FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=gdpoore@eng.ucsd.edu   # Where to send mail	
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --mem=16gb                   # Job Memory
#SBATCH --time=00:45:00             # Time limit hrs:min:sec
#SBATCH --output=/home/gdpoore/cluster/logs/array_%A_%a.log    # Standard output and error log
#SBATCH -a 0-999%25                # Array range
### RUN LIKE THIS: sbatch array_grep_fungal_OGUs.sh "file_list_alignments_qiita_rep200_wgs_and_rna_4.txt"

filelist=$1

in_dir=/projects/tcga-data/wgs_extract/alignments_qiita_rep200_wgs_and_rna
out_dir=/projects/tcga-data/mycobiome_coverage_analysis/old_extracted_sam_files_fungi_all_qiita_samples
filelist_dir=/projects/tcga-data/mycobiome_coverage_analysis/file_lists
fungi_old_gOTU_file=/projects/tcga-data/rna_extract/norm_fungi_rep200_gID_only.txt
cd ${in_dir}

readarray -t filenameArray < ${filelist_dir}/${filelist}
filename=$(echo ${filenameArray[${SLURM_ARRAY_TASK_ID}]} | tr -d '\r')
no_root=${filename##*/}
base_name=${no_root%.sam.xz*}

source activate ebi_sra_importer

unxz -c ${in_dir}/${base_name}.sam.xz | grep -f ${fungi_old_gOTU_file} > ${out_dir}/filt_fungi_old_${base_name}.sam