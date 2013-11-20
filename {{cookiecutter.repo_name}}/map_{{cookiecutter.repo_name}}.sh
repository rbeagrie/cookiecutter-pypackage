#!/bin/bash
#$ -l h_vmem=12G
#$ -cwd
set -e

echo "TopHat version:"
{{ cookiecutter.tophat_binary }} --version

READS="`./get_reads.sh read_symlinks/*_R1*` `./get_reads.sh read_symlinks/*_R2*`"
echo 
python check_reads.py $READS
echo "Finished checking reads, starting the mapping."

{{ cookiecutter.tophat_binary }} -p {{ cookiecutter.num_processes }} -r {{ cookiecutter.mate_inner_dist }} --mate-std-dev {{ cookiecutter.mate_std_dev }} --library-type {{ cookiecutter.library_type }} -G {{ cookiecutter.gtf_file }} -o mapping_results {{ cookiecutter.bowtie_index }} $READS

