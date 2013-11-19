#!/bin/bash
#$ -l h_vmem=12G
#$ -cwd

READS=""

{{ cookiecutter.tophat_binary }} -p {{ cookiecutter.num_processes }} -r {{ cookiecutter.mate_inner_dist }} --mate-std-dev {{ cookiecutter.mate_std_dev }} --library-type {{ cookiecutter.library_type }} -G {{ cookiecutter.gtf_file }} -o mapping_result {{ cookiecutter.bowtie_index }} $READS

