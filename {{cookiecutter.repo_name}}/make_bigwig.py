#! /usr/bin/env python

import sys
import argparse
import os
from doit.cmd_base import ModuleTaskLoader
from doit.doit_cmd import DoitMain


parser = argparse.ArgumentParser(description='Create BigWigs from RNA-seq data.')
parser.add_argument('input_bamfiles', metavar='INPUT_BAMFILE', nargs='+', help='Input bam file(s) to create BigWigs from.')
parser.add_argument('-g','--genome_file', metavar='GENOME_FILE', required=True, help='File containing chromosome names and lengths')
parser.add_argument('-s','--stranded',action='store_true',
                    help='Create two BigWigs, one for the plus and one for the minus strand (default: one file for both strands)')

def with_extension(input_file, extension): 
    output_dir = os.path.dirname(input_file)
    basename = os.path.splitext(os.path.basename(input_file))[0]
    return os.path.join(output_dir, basename + extension)


def sort_bam():
    """Sort a bamfile by genomic position"""

    for input_bamfile in args.input_bamfiles:
        yield {
            'basename' : 'Sorting input bam',
            'name'     : input_bamfile,
            'actions'  : ['samtools sort %(dependencies)s %(targets)s',
                          'mv %(targets)s.bam %(targets)s'],
            'targets'  : [with_extension(input_bamfile, ".sorted.bam")],
            'file_dep' : [input_bamfile],
              }

def get_plus_bamfile():
    """Extract reads mapping to the positive strand of a bamfile"""

    for input_bamfile in args.input_bamfiles:
        yield {
            'basename' : 'Extracting plus strand reads.',
            'name'     : input_bamfile,
            'actions'  : ['bash -i -c "samtools merge %(targets)s <(samtools view -hb -f 128 -F 16 %(dependencies)s) <(samtools view -hb -f 80 %(dependencies)s)"'],
            'targets'  : [with_extension(input_bamfile, ".plus.bam")],
            'file_dep' : [with_extension(input_bamfile, ".sorted.bam")],
        }

def get_minus_bamfile():
    """Extract reads mapping to the negative strand of a bamfile"""

    for input_bamfile in args.input_bamfiles:
        yield {
            'basename' : 'Extracting minus strand reads',
            'name'     : input_bamfile,
            'actions'  : ['bash -i -c "samtools merge %(targets)s <(samtools view -hb -F 16 -f 64 %(dependencies)s) <(samtools view -hb -f 144 %(dependencies)s)"'],
            'targets'  : [with_extension(input_bamfile, ".minus.bam")],
            'file_dep' : [with_extension(input_bamfile, ".sorted.bam")],
        }

def get_plus_bedgraph():
    """Convert positive strand bamfiles to bedgraphs"""

    for input_bamfile in args.input_bamfiles:
        yield {
            'basename' : 'Making bedGraph for plus strand',
            'name'     : input_bamfile,
            'actions'  : ['genomeCoverageBed -bg -split -ibam %(dependencies)s -g %(genome_file)s > %(targets)s'],
            'targets'  : [with_extension(input_bamfile, ".plus.bedgraph")],
            'file_dep' : [with_extension(input_bamfile, ".plus.bam")],
        }

def get_plus_bigwig():
    """Convert positive strand bedGraphs to BigWigs"""

    for input_bamfile in args.input_bamfiles:
        yield {
            'basename' : 'Converting plus strand bedGraph to BigWig',
            'name'    : input_bamfile,
            'actions' : ['bedGraphToBigWig %(dependencies)s %(genome_file)s %(targets)s'],
            'targets' : [with_extension(input_bamfile, ".plus.bw")],
            'file_dep' : [with_extension(input_bamfile, ".plus.bedgraph")],
        }

def get_minus_bedgraph():
    """Convert minus strand bamfiles to bedGraphs"""

    for input_bamfile in args.input_bamfiles:
        yield {
            'basename' : 'Making bedGraph for minus strand',
            'name'     : input_bamfile,
            'actions'  : ['genomeCoverageBed -bg -split -ibam %(dependencies)s -g %(genome_file)s > %(targets)s.temp',
                         'awk \'{print $1"\t"$2"\t"$3"\t-"$4}\' %(targets)s.temp > %(targets)s',
                         'rm %(targets)s.temp'],
            'targets'  : [with_extension(input_bamfile, ".minus.bedgraph")],
            'file_dep' : [with_extension(input_bamfile, ".minus.bam")],
        }

def get_minus_bigwig():
    """Convert minus strand bedGraphs to BigWigs"""

    for input_bamfile in args.input_bamfiles:
        yield {
            'basename' : 'Converting minus strand bedGraph to BigWig',
            'name'     : input_bamfile,
            'actions'  : ['bedGraphToBigWig %(dependencies)s %(genome_file)s %(targets)s'],
            'targets'  : [with_extension(input_bamfile, ".minus.bw")],
            'file_dep' : [with_extension(input_bamfile, ".minus.bedgraph")],
        }

def get_unstranded_bedgraph():
    """Convert unstranded bamfiles to bedGraphs"""

    for input_bamfile in args.input_bamfiles:
        yield {
            'basename' : 'Converting input bam file to bedGraph',
            'name'     : input_bamfile,
            'actions'  : ['genomeCoverageBed -bg -split -ibam %(dependencies)s -g %(genome_file)s > %(targets)s'],
            'targets'  : [with_extension(input_bamfile, ".bedgraph")],
            'file_dep' : ['input_bamfile'],
        }

def get_unstranded_bigwig():
    """Convert unstranded bamfiles to bedGraphs"""

    for input_bamfile in args.input_bamfiles:
        yield {
            'basename' : 'Converting bedGraph file to BigWig',
            'name'     : input_bamfile,
            'actions'  : ['bedGraphToBigWig %(dependencies)s %(genome_file)s %(targets)s'],
            'targets'  : [with_extension(input_bamfile, ".bw")],
            'file_dep' : [with_extension(input_bamfile, ".bedgraph")],
        }

def main(args):

    tasks_to_run = { 'task_sort_bam' : sort_bam }

    stranded_tasks = { 'task_plus_bamfile'   : get_plus_bamfile,
                       'task_minus_bamfile'  : get_minus_bamfile,
                       'task_plus_bedgraph'  : get_plus_bedgraph,
                       'task_minus_bedgraph' : get_minus_bedgraph,
                       'task_plus_bigwig'    : get_plus_bigwig,
                       'task_minus_bigwig'   : get_minus_bigwig }

    unstranded_tasks = { 'task_bedgraph' : get_unstranded_bedgraph,
                         'task_bigwig'   : get_unstranded_bigwig }

    if args.stranded:
        tasks_to_run.update(stranded_tasks)
    else:
        tasks_to_run.update(unstranded_tasks)
    
    sys.exit(DoitMain(ModuleTaskLoader(tasks_to_run)).run([]))

if __name__ == "__main__":

    args = parser.parse_args()

    main(args)
