# Pan-cancer mycobiome atlas

This is the GitHub repository for all TCGA-centric analyses of the pan-cancer mycobiome atlas paper by Narunsky-Haziza, Sepich-Poore, Livyatan _et al.,_ 2022. Data and code provided herein cover the fungal and bacterial analysis of 17,401 patient tissue, blood, and plasma samples across 35 cancer types in four independent cohorts. 

**Please note** that this repository is divided into multiple subdirectories, each with their respective files and/or instructions, as shown below.

## 1. Dockerized host depletion pipeline (Under `/Docker_host_depletion_pipeline`)

The following files are listed in this subdirectory:
1. `Dockerfile` —> This text file is used to build the Docker container on your computer (instructions are below)
2. `ebi_sra_importer.yml` —> This yml file contains the conda environment and its dependences that are run in the Docker container
3. `host_deplete_script.sh` —> This bash script is used to run the host depletion within the Docker container on any new file. It accepts 3 arguments (in this order): (a) Original bam file to host deplete, (b) Minimap2 human database in .mmi format, (c) Integer number of cpus you would like to use on your server or machine.

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

If you use data or software from this repository, please cite the following two papers:

- Narunsky-Haziza, Sepich-Poore, Livyatan _et al._
```
{Citation for this paper forthcoming}
```
As well as Poore, Kopylova _et al._:
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
