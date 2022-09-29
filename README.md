# Pan-cancer mycobiome atlas

This is the GitHub repository for all TCGA, Hopkins, and UCSD analyses of the pan-cancer mycobiome atlas paper by Narunsky-Haziza, Sepich-Poore, Livyatan _et al.,_ 2022. Data and code provided herein cover the fungal and bacterial analysis of 17,401 patient tissue, blood, and plasma samples across 35 cancer types in four independent cohorts.

<p align="center">
<img src="Figures/Graphical_abstract/Mycobiome-graphical-abstract-final.png" width="700">
</p>

**Please note** that this repository is divided into multiple subdirectories. They are detailed below in the `Organization of files` section.

## Related repositories

There are several other repositories related to this paper, including the following:

- The ITS2 pipeline for amplicon-based analyses: https://github.com/microbiofunc/ITS2-pipeline
- Fungal phylogenomics of reconstructed bins: https://github.com/stajichlab/Tumor_mycobiome_2022
- Aggregated and per-sample genome coverage calculations: https://github.com/ucsd-cmi/zebra_filter

## Installation

On the top right of the main repository page, click the `Code ▼` button followed by the `Download ZIP` button to download the entire repository.

Alternatively, or if using `R`, you can download the repository using:
```
library(devtools)
install_github("gregpoore/mycobiome")
```

## Organization of files

### `R` files, scripts, and subfolders:

The `R` files are enumerated and provided approximately in the same order as TCGA analyses detailed in the paper. The following subdirectories contain important information for running the `R` analyses:

- All of the input data used for the 17 main R scripts in the base directory are provided in the `/Input_data` subfolder. 
- Supporting data files (e.g., taxonomy tables) are provided under the `/Supporting_data` subfolder.
- Data generated by earlier scripts that are used by later scripts are saved under `/Interm_data`. 
- Scripts that were run on a Slurm compute cluster using processed versions of the data are located under `/Supporting_scripts`.
- Figures that were generated while processing the scripts were saved under `/Figures` (note the multiple subdirectories therein).
- In certain cases (e.g., ANCOM-BC differential abundance tables), data underlying certain Figures were saved under `/Figures_data`.

The 17 main `R` scripts in the base directory have the following brief descriptions:

- `00-Functions.R` → This script contains all major functions used throughout the rest of the `R` scripts
- `01R-Merge-WGS-RNA-data-decontaminate-batch-correct.R` → This script merged the TCGA WGS and RNA-Seq datasets (based on their [Qiita](https://qiita.ucsd.edu/) outputs), followed by data decontamination and batch correction.
- `02R-Calculate-fungi-vs-bacteria-read-distributions.R` → This script calculated the relative abundances of fungi vs. bacteria with and without genome size correction, and the fungal vs. bacterial read percentages.
- `03R-Intersect-TCGA-data-with-Weizmann-data.R` → This script intersects the TCGA and WIS fungal and bacterial data at every taxa level.
- `04R-Prepare-TCGA-data-for-Qiime-and-plot-alpha-diversity.R` → This script saves data for [Qiime 2](https://qiime2.org/) alpha and beta diversity analyses of TCGA fungal data. **Note:** The data saved here were processed using `/Qiime_data_and_scripts_resubmission_version/qiime2_mycobiome_tcga_analyses_resubmission_version.ipynb`.
- `05R-Prepare-TCGA-data-for-machine-learning.R` → This script does a lot of data formatting for most of the machine learning analyses presented in the paper. Typically, the formatted data are saved as `.RData` files, which are then called by Slurm-based scripts (stored under `/Supporting_scripts`), and the results of those Slurm-based runs are then plotted in this script.
- `06R-Perform-machine-learning-on-Weizmann-data.R` → This script formats the WIS data for machine learning and does a preliminary version of it (although the main version of the machine learning is under `/Supporting_scripts/S16R-ML-fungi-10k-rep1-weizmann.R`.
- `07R-TCGA-compare-tumor-vs-normal-bray-curtis.R` → This script calculates a Bray-Curtis-based PCoA for tumor vs. NAT samples using cancer types that overlapped with the WIS cohort. The resulting output comprises a 2D and 3D plot.
- `08R-TCGA-alpha-diversity-correlation-fungi-vs-bacteria.R` → This script performs rarefaction on the entire dataset and calculates correlations among fungal and bacterial diversities across all cancer types.
- `09R-Prepare-TCGA-data-for-MMvec-and-ANCOM-BC.R` → This script formats the data for running MMvec to identify fungal-bacterial-immune co-occurrences (see `/MMvec-cooccurrence-analyses-resubmission` for more) and runs ANCOM-BC differential abundance on fungi and bacteria in TCGA.
- `10R-Control-validation-analyses-TCGA.R` → This script performs control analyses, including stratified splits in TCGA with cross-testing and permuted machine learning data analyses.
- `11R-Plasma-validation-cohort-UCSD.R` → This script re-analyses plasma microbiome data generated by [Poore and Kopylova _et al._ 2020 Nature](https://www.nature.com/articles/s41586-020-2095-1) for fungal analytes.
- `12R-Plasma-validation-cohort-Hopkins.R` → This script re-analyses shallow WGS multi-cancer data generated by [Cristiano _et al._ 2020 Nature](https://www.nature.com/articles/s41586-019-1272-6) for fungal and other microbial analytes.
- `13R-ML-fungi-vs-bacteria.R` → This script examines the synergy of using fungal and bacterial biomarkers in combination for pan-cancer machine learning.
- `14R-ML-fungi-cancer-stage.R` → This script performs machine learning to discriminate early vs. late stage cancers.
- `15R-ML-WIS-samples.R` → This script performs the remaining machine learning analyses on WIS data that were not covered in `06R-Perform-machine-learning-on-Weizmann-data.R`.
- `16R-Addressing-remaining-reviewer-comments.R` → This script performs several analyses to respond to reviewer requests.

### Qiime alpha and beta diversity analyses:

The Qiime-based analyses, including input data and script, are found under the `/Qiime_data_and_scripts_resubmission_version` subfolder. The `CLI` commands run are listed therein under `qiime2_mycobiome_tcga_analyses_resubmission_version.ipynb`.

### MMvec fungal-bacterial-immune co-occurrence analyses:

The MMvec-based analyses, including input data and scripts, are found under the `/MMvec-cooccurrence-analyses-resubmission` subfolder. A separate `README.md` file is listed within that subdirectory explaining its contents.

### Other directories:

- `/EukDetect` → This directory contains code and resultant tables for rerunning TCGA using [EukDetect](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01015-y).
- `/Coverage_analysis` → This directory contains several scripts involved in calculating fungal genome coverages.

## Docker host depletion pipeline (`/Docker_host_depletion_pipeline`)

For these analyses, we re-aligned all non-human reads against a uniform reference genome (GRCh38+PhiX). Due to the large amount of data being reprocessed, the host-depletion pipeline was optimized for speed. For others to use, we have packaged the host depletion steps into a Docker container, described below.

The following files are listed in this subdirectory:
1. `Dockerfile` —> This text file is used to build the Docker container on your computer (instructions are below)
2. `ebi_sra_importer.yml` —> This yml file contains the conda environment and its dependences that are run in the Docker container
3. `host_deplete_script.sh` —> This bash script is used to run the host depletion within the Docker container on any new file. It accepts 3 arguments (in this order): (a) Original bam file to host deplete, (b) Minimap2 human database in .mmi format, (c) Integer number of cpus you would like to use on your server or machine.

**Note:** The minimap2 database (.mmi file) is located on figshare here: https://doi.org/10.6084/m9.figshare.19719418.v1

### Building the Docker container:

You can build the Docker container with the following command (please substitute the `<…>` text with your own text; note the `.` at the end of the command as well): `docker build --no-cache -t <MY_CONTAINER_NAME> .`

It will take some time (~20-30 minutes) to build the container and occupy several gigabytes of space. Once the container has been built, you can interactively run it with the following command: `docker run -ti <MY_CONTAINER_NAME>`

Once the container is running, you can use the “host_deplete_script.sh” script to run the host depletion on any bam file you would like. It will then save the host depleted file as two fastq files with the following suffixes: `<BASE_NAME>.R1.trimmed.fastq.gz`, `<BASE_NAME>.R2.trimmed.fastq.gz`

The script within the Docker container is running the following one-liner piped command:
```
samtools view -f 4 -O BAM $in_dir/$filename |
samtools bam2fq - |
fastp -l 45 --stdin -w $cpus --stdout --interleaved_in |
minimap2 -ax sr -t $cpus $db - |
samtools fastq -@ $cpus -f 12 -F 256 - -1 $out_dir/$base_name.R1.trimmed.fastq.gz -2 $out_dir/$base_name.R2.trimmed.fastq.gz
```
where `$cpus` and `$db` denote the number of compute cores and a precomputed Minimap2 reference database (as a .mmi file), respectively.

**Note 1:** The minimap2 database for this study included the [GRCh38.p7 human genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.33/) and the Phi X 174 viral genome.

**Note 2:** The host depletion automatically removes any reads less than 45 bp for quality control. If you would like to change this, you can alter the `-l 45` argument to fastp on line 9 of the `host_deplete_script.sh`.

**Note 3:** If you prefer to avoid using a Docker container, you can create a new conda environment using the `ebi_sra_importer.yml` file on any machine or server with Conda installed. You can then adapt the `host_deplete_script.sh` script as you like to run within the Conda environment on any bam file that you would like to host deplete.

After the host depletion is complete: You should have paired R1 and R2 fastq files for each host depleted sample, which you can then upload to [Qiita](https://qiita.ucsd.edu/) to perform taxonomy calling using [Woltka](https://github.com/qiyunzhu/woltka).


## Citation

If you use data or software from this repository, please cite the following <ins>two</ins> papers:

- Narunsky-Haziza, Sepich-Poore, Livyatan _et al._
```
{Citation for this paper forthcoming}
```
- [Poore, Kopylova _et al._ 2020. _Nature_](https://www.nature.com/articles/s41586-020-2095-1):
```
@article{poore2020microbiome,
  title={Microbiome analyses of blood and tissues suggest cancer diagnostic approach},
  author={Poore, Gregory D and Kopylova, Evguenia and Zhu, Qiyun and Carpenter, Carolina and Fraraccio, Serena and Wandro, Stephen and Kosciolek, Tomasz and Janssen, Stefan and Metcalf, Jessica and Song, Se Jin and others},
  journal={Nature},
  volume={579},
  number={7800},
  pages={567--574},
  year={2020},
  publisher={Nature Publishing Group}
}
```


## License
**Please carefully note the license terms surrounding commercial-use**

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
