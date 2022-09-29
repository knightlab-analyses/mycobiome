
## Download final TCGA mycobiome files
The final raw and batch-corrected TCGA files that were used for machine learning are located in this. This includes fungal data for 14,495 samples. The files are the following:

- A taxonomy table for the rep200 database used throughout out paper (`taxonomy_table_rep200.tsv`).
- A metadata table for the 14,495 TCGA samples (`metadata_fungi_14495samples.tsv`). *Note* that the TCGA case IDs are under the `tcga_case_id` column.
- A raw count table for 224 decontaminated fungi (`count_data_fungi_decontaminated_raw.tsv`).
- A batch-corrected (using Voom-SNM) table for 224 decontaminated fungi (`count_data_fungi_decontaminated_voom_snm_corrected.tsv`).
