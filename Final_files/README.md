
## Download final TCGA decontaminated fungi files
The final raw and batch-corrected TCGA files that were used for machine learning are located in this. This includes fungal data for 14,495 samples. The files are the following:

- A taxonomy table for the rep200 database used throughout out paper (`taxonomy_table_rep200.tsv`).
- A metadata table for the 14,495 TCGA samples (`metadata_fungi_14495samples.tsv`). *Note* that the TCGA case IDs are under the `tcga_case_id` column.
- A raw count table for 224 decontaminated fungi (`count_data_fungi_decontaminated_raw.tsv`).
- A batch-corrected (using Voom-SNM) table for 224 decontaminated fungi (`count_data_fungi_decontaminated_voom_snm_corrected.tsv`).

## Download TCGA fungi+bacterial data
Since the release of the paper, we have received requests to release the bacterial data alongside the fungi information. For now, we are conservatively posting abundances of Weizmann-overlapping bacteria (cf. `Nejman et al. 2020. Science`) and fungi (this paper, `Narunsky-Haziza et al. 2022. Cell`) at the genus and species levels, since these microbes have multiple layers of evidence in independent cohorts. These count files incorporate the updated host-depletion and read QC methods described in this paper. The following files contain Weizmann-overlapping bacteria and fungi in TCGA:

- TCGA raw genus-level abundances for WIS-overlapping bacteria and fungi (14,494 non-zero samples) (`count_data_genus_raw_WIS_overlapping_fungi_bacteria_14494samples.tsv`)
- Matching TCGA metadata for genus-level data (14,494 non-zero samples) (`metadata_genus_WIS_overlapping_fungi_bacteria_14494samples.tsv`)
- TCGA raw species-level abundances for WIS-overlapping bacteria and fungi (12,773 non-zero samples) (`count_data_species_raw_WIS_overlapping_fungi_bacteria_12773samples.tsv`)
- Matching TCGA metadata for species-level data (12,773 non-zero samples) (`metadata_species_WIS_overlapping_fungi_bacteria_12773samples.tsv`)
- Matching taxonomy table at species-level for WIS-overlapping bacteria and fungi (`taxonomy_table_WIS_overlapping_fungi_bacteria.tsv`)

**Note #1:** Subsetting TCGA samples to WIS-overlapping microbes leads to some sample dropout, which is why the sample counts in these files are slightly lower than the 14,495 number posted at the top. 

**Note #2:** The raw counts/abundances posted above do _not_ correct for batch effects in TCGA. One should either mitigate the batch effects using batch correction (e.g., Voom-SNM, ComBat-Seq, etc.) or subset the samples to groups not affected by individual batches (i.e., a single sequencing center, experimental strategy [WGS or RNA-Seq], and sequencing platform).

**Note #3:** If using any of these files, please cite *both* `Nejman et al. 2020. Science` and `Narunsky-Haziza et al. 2022. Cell` (this paper).

**Note #4:** Completely raw data, including the host-depleted fastq files, can be accessed at the respective Qiita project study IDs: 13722 (TCGA WGS), 13767 (TCGA RNA-Seq). 
