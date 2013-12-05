#!/bin/bash
#$ -l h_vmem=12G
#$ -cwd
#$ -V
set -e

echo "TopHat version:"
{{ cookiecutter.tophat_binary }} --version

READS="`./get_reads.sh read_symlinks/*_R1*` `./get_reads.sh read_symlinks/*_R2*`"
echo 
python check_reads.py $READS
echo "Finished checking reads, starting the mapping."

{{ cookiecutter.tophat_binary }} -p {{ cookiecutter.num_processes }} -r {{ cookiecutter.mate_inner_dist }} --mate-std-dev {{ cookiecutter.mate_std_dev }} --library-type {{ cookiecutter.library_type }} -G {{ cookiecutter.gtf_file }} -o mapping_results {{ cookiecutter.bowtie_index }} $READS

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
chmod +x fetchChromSizes
./fetchChromSizes {{ cookiecutter.genome_db }} > {{ cookiecutter.genome_db }}.chrom.sizes

{{ cookiecutter.python_path }} -u make_bigwig.py -g {{ cookiecutter.genome_db }}.chrom.sizes {{ cookiecutter.stranded }} mapping_results/accepted_hits.bam
