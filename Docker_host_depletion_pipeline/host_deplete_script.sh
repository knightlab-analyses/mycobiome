#!/bin/bash -l
filename=$1
db=$2
cpus=$3
no_root=${filename##*/}
base_name=${no_root%.bam*}
samtools view -f 4 -O BAM ${filename} |
samtools bam2fq - |
fastp -l 45 --stdin -w ${cpus} --stdout --interleaved_in |
minimap2 -ax sr -t ${cpus} ${db} - |
samtools fastq -@ ${cpus} -f 12 -F 256 - -1 ${base_name}.R1.trimmed.fastq.gz -2 ${base_name}.R2.trimmed.fastq.gz