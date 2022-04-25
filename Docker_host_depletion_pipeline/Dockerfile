FROM continuumio/miniconda:latest
MAINTAINER greg poore
ADD ebi_sra_importer.yml /code/ebi_sra_importer.yml
ADD host_deplete_script.sh /host_deplete_script.sh
RUN conda env create -f /code/ebi_sra_importer.yml
RUN echo "conda activate $(head -1 /code/ebi_sra_importer.yml | cut -d' ' -f2)" >> ~/.bashrc
ENV PATH /opt/conda/envs/$(head -1 /code/ebi_sra_importer.yml | cut -d' ' -f2)/bin:$PATH
ENV CONDA_DEFAULT_ENV $(head -1 /code/ebi_sra_importer.yml | cut -d' ' -f2)
ENV PATH /opt/conda/envs/ebi_sra_importer/bin:/opt/conda/condabin:/opt/conda/envs/$(head -1 /code/ebi_sra_importer.yml | cut -d -f2)/bin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/conda/envs/ebi_sra_importer/bin/samtools:/opt/conda/envs/ebi_sra_importer/bin/fastp:/opt/conda/envs/ebi_sra_importer/bin/minimap