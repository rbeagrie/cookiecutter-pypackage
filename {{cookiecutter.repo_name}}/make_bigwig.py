#! /usr/bin/env python

from doit.task import dict_to_task
from doit.cmd_base import TaskLoader
from doit.doit_cmd import DoitMain
import argparse
import os
import sys

parser = argparse.ArgumentParser(description='Create BigWigs from RNA-seq data.')
parser.add_argument('input_bamfile', metavar='INPUT_BAMFILE', help='Input bam file to create BigWigs from.')
parser.add_argument('-g','--genome_file', metavar='GENOME_FILE', help='File containing chromosome names and lengths')
parser.add_argument('-s','--stranded',action='store_true',
                    help='Create two BigWigs, one for the plus and one for the minus strand (default: one file for both strands)')

args = parser.parse_args()

def get_doit_tasks(input_bamfile):
    """Take the path to one bamfile and return the doit tasks which would convert it to a BigWig"""

    output_dir = os.path.dirname(input_bamfile)
    basename = os.path.splitext(os.path.basename(input_bamfile))[0]

    def with_extension(extension): return os.path.join(output_dir, basename + extension)

    tasks_to_run = []

    tasks_to_run.append({
                'name' : 'Sorting input bam',
                'actions' : ['samtools sort %(dependencies)s %(targets)s',
                          'mv %(targets)s.bam %(targets)s'],
                'targets' : [with_extension(".sorted.bam")],
                'file_dep' : [input_bamfile],
    })

    if args.stranded:
        tasks_to_run.append({
                    'name' : 'Extracting plus strand reads.',
                    'actions' : ['bash -i -c "samtools merge %(targets)s <(samtools view -hb -f 128 -F 16 %(dependencies)s) <(samtools view -hb -f 80 %(dependencies)s)"'],
                    'targets' : [with_extension(".plus.bam")],
                    'file_dep' : [with_extension(".sorted.bam")],
        })

        tasks_to_run.append({
                    'name' : 'Extracting minus strand reads',
                    'actions' : ['bash -i -c "samtools merge %(targets)s <(samtools view -hb -F 16 -f 64 %(dependencies)s) <(samtools view -hb -f 144 %(dependencies)s)"'],
                    'targets' : [with_extension(".minus.bam")],
                    'file_dep' : [with_extension(".sorted.bam")],
        })

        tasks_to_run.append({
                    'name' : 'Making bedGraph for plus strand',
                    'actions' : ['genomeCoverageBed -bg -split -ibam %(dependencies)s -g %(genome_file)s > %(targets)s'],
                    'targets' : [with_extension(".plus.bedgraph")],
                    'file_dep' : [with_extension(".plus.bam")],
        })

        tasks_to_run.append({
                    'name' : 'Converting plus strand bedGraph to BigWig',
                    'actions' : ['bedGraphToBigWig %(dependencies)s %(genome_file)s %(targets)s'],
                    'targets' : [with_extension(".plus.bw")],
                    'file_dep' : [with_extension(".plus.bedgraph")],
        })

        tasks_to_run.append({
                    'name' : 'Making bedGraph for minus strand',
                    'actions' : ['genomeCoverageBed -bg -split -ibam %(dependencies)s -g %(genome_file)s > %(targets)s'],
                    'targets' : [with_extension(".minus.bedgraph")],
                    'file_dep' : [with_extension(".minus.bam")],
        })

        tasks_to_run.append({
                    'name' : 'Converting minus strand bedGraph to BigWig',
                    'actions' : ['bedGraphToBigWig %(dependencies)s %(genome_file)s %(targets)s'],
                    'targets' : [with_extension(".minus.bw")],
                    'file_dep' : [with_extension(".minus.bedgraph")],
        })

    else:
        tasks_to_run.append({
                    'name' : 'Converting input bam file to bedGraph',
                    'actions' : ['genomeCoverageBed -bg -split -ibam %(dependencies)s -g %(genome_file)s > %(targets)s'],
                    'targets' : [with_extension(".bedgraph")],
                    'file_dep' : ['input_bamfile'],
        })

        tasks_to_run.append({
                    'name' : 'Converting bedGraph file to BigWig',
                    'actions' : ['bedGraphToBigWig %(dependencies)s %(genome_file)s %(targets)s'],
                    'targets' : [with_extension(".bw")],
                    'file_dep' : [with_extension(".bedgraph")],
        })

    return tasks_to_run

class MyTaskLoader(TaskLoader):
    def __init__(self, tasks_to_run, args):
        self.tasks_to_run = tasks_to_run
        self.doit_db = os.path.join(self.get_db_dir(args),'.doit.db')
        self.substitute_tasks(args)
        super(TaskLoader, self).__init__()

    def get_db_dir(self, args):

        return os.path.dirname(args.input_bamfile)

    def substitute_tasks(self, args):

        subs_dictionary = vars(args)
        subs_dictionary.update({'dependencies' : '%(dependencies)s',
                                'targets' : '%(targets)s'})
        print subs_dictionary

        for task in self.tasks_to_run:
            task['actions'] = [ action % subs_dictionary for action in task['actions'] ]

    #@staticmethod
    def load_tasks(self, cmd, opt_values, pos_args):
        task_list = [ dict_to_task(my_task) for my_task in self.tasks_to_run ]
        config = {'verbosity': 2,
                  'dep_file' : self.doit_db
                 }
        return task_list, config

def main(args):

    tasks_to_run = get_doit_tasks(args.input_bamfile)

    sys.exit(DoitMain(MyTaskLoader(tasks_to_run, args)).run([]))

if __name__ == "__main__":

    args = parser.parse_args()

    main(args)
