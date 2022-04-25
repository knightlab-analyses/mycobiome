# => processing via fastp
# grep -v '#' /projects/tcga-data/EukDetect/prep_artifact_ids.txt | while IFS=' ' read pid aid ; do mkdir -p /panfs/panfs1.ucsd.edu/panscratch/antoniog/tcga/${pid}-${aid}/ ; for f in `ls /qmounts/qiita_data/per_sample_FASTQ/${aid}/*.R1.*fastq.gz`; do fname=`basename ${f}`; echo fastp -l 75 -w 15 -i ${f} -I ${f/.R1./.R2.} -o /panfs/panfs1.ucsd.edu/panscratch/antoniog/tcga/${pid}-${aid}/${fname} -O /panfs/panfs1.ucsd.edu/panscratch/antoniog/tcga/${pid}-${aid}/${fname/.R1./.R2.} ; done ; done > /projects/tcga-data/EukDetect/fastp_commands.sh
# conda create -n fastp -c bioconda fastp
# echo "source ~/.bash_profile; conda activate fastp; sh /projects/tcga-data/EukDetect/fastp_commands.sh" | qsub -N fastp -l mem=20gb -l walltime=300:00:00 -l nodes=1:ppn=15

# => run from /panfs/panfs1.ucsd.edu/panscratch/antoniog/tcga
#    move empty files
#    check file counts
# for d in `ls -d 10*`; do count=`ls -s ${d}/*gz | wc -l`; echo ${d} $count;  done
# for d in `ls -d 10*`; do echo ${d}; ls -s ${d}/*gz | while IFS=' ' read z f ; do if [ "$z" -le "48" ]; then mkdir -p ${d}/empty/; mv $f ${d}/empty/; fi;  done ; done

# => run EukDetect:
#    creating config files
# for d in `ls -d 10*`; do outdir=${PWD}/Euk-${d}; echo $d; rfile=`ls $d/*.gz | head -n 1`; readlen=`gzip -dc ${rfile} | head -n 10000 | awk '{ if (NR%4==2){count++; bases += length}} END{printf "%3.0f\n", bases/count}'`; mkdir -p ${outdir}; main_f=`cat /projects/tcga-data/EukDetect/base_configfile.yml`; files=`ls ${d}/*.R1.trimmed.fastq.gz | awk  -F '/' '{sub(/.R1.trimmed.fastq.gz/, ""); print "  "$NF":"}'`; echo -e "fq_dir: \"${PWD}/${d}\"\noutput_dir: \"${outdir}\"\nreadlen: ${readlen}\n${main_f}\n\nsamples:\n${files}" > ${outdir}/configfile.yml ; done

#    submitting jobs
# for f in `ls */configfile*`; do dname=`dirname $f`; echo "source ~/.bash_profile; conda activate eukdetect; cd $PWD/${dname}; eukdetect --mode runall --configfile $PWD/${f} --cores 16" | qsub -l walltime=100:00:00 -l mem=100g -l nodes=1:ppn=16 -N $dname; done

# ** There is one folder that has way to mixed sequence length so we need to fix it
# cd 10653-115226
# for f in `ls *.R1.*`; do readlen=`gzip -dc ${f} | head -n 10000 | awk '{ if (NR%4==2){count++; bases += length}} END{printf "%3.0f\n", bases/count}'`; echo $readlen $f; done | sort -n > ../Euk-10653-115226/len_count.txt
# cd ../Euk-10653-115226
# body=`head -n 12 configfile.yml`; files=`head -n 236 len_count.txt | awk '{sub(/.R1.trimmed.fastq.gz/, ""); print "  "$NF":"}'`; echo -e "${body/: 150/: 100}\n${files}" > configfile_100.yml
# files=`tail -n 14 len_count.txt | awk '{sub(/.R1.trimmed.fastq.gz/, ""); print "  "$NF":"}'`; echo -e "${body}\n${files}" > configfile_150.yml
# echo "source ~/.bash_profile; conda activate eukdetect; cd $PWD; eukdetect --mode runall --configfile $PWD/configfile_100.yml --cores 16" | qsub -l walltime=100:00:00 -l mem=100g -l nodes=1:ppn=16 -N Euk100;
# echo "source ~/.bash_profile; conda activate eukdetect; cd $PWD; eukdetect --mode runall --configfile $PWD/configfile_150.yml --cores 16" | qsub -l walltime=100:00:00 -l mem=100g -l nodes=1:ppn=16 -N Euk150;

# => generating tables!!

# python /projects/tcga-data/EukDetect/generate_tables.py
# for f in `ls tables/*.tsv | grep -v taxonomy`; do echo $f; biom convert -i $f -o ${f/tsv/biom} --table-type="OTU table" --to-hdf5; qiime tools import --input-path ${f/tsv/biom} --type 'FeatureTable[Frequency]' --output-path ${f/tsv/qza}; done
