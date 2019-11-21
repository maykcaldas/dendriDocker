#!/usr/bin/python3
#!-*- coding: utf8 -*-

'''
Software general documentation

gromacsBuilding was made using python 3.5.2

'''

from writeMDP import *
from writeRun import *
import argparse
import time
import workflow

def write_submission(file_name, job_name, job_file):
    if os.path.isfile(file_name):
        print('!!!Backing up the existing job file!!!')
        os.rename(file_name, 'bck.'+file_name)
    job=open(file_name,'w')
 
    job.write('#!/bin/bash\n')
    job.write('#PBS -l select=1:ncpus=48:mpiprocs=8\n')
    job.write('#PBS -l walltime=999:00:00\n')
    job.write('#PBS -j oe\n')
    job.write('#PBS -M maykcaldas@gmail.com\n')
    job.write('#PBS -m bea\n')
    job.write('#PBS -V\n')
    job.write('#PBS -N {job_name}\n'.format(job_name=job_name))
    job.write('\n')
    job.write('# load modules\n')
    job.write('module load openmpi-gnu/2.1.1\n')
    job.write('source /scratch/60061a/plumed2/sourceme.sh\n')
    job.write('# change directory\n')
    job.write('cd $\{PBS_O_WORKDIR\}\n')
    job.write('# environment (if necessary)\n')
    job.write('export PLUMED_NUM_THREADS=1\n')
    job.write('export OMP_NUM_THREADS=6\n')
    job.write('\n')
    job.write('# run\n')
    job.write('{job_file}'.format(job_file=job_file))


def write_log(log_file, workFlow):

    print("The workFlow is:")
    flow=""
    for step in workFlow:
        flow+="-----> {0} ".format(step[0])
    print(flow)

    if os.path.isfile(log_file):
        print('!!!Backing up the existing log file!!!')
        os.rename(log_file, 'bck.'+log_file)
    log=open(log_file,'a')

    for step in workFlow:
        log.write(step[0])
        log.write("\t mdp input\n")
        log.write("\t {0}\n".format(step[1]))
        log.write("\t run input\n")
        log.write("\t {0}\n".format(step[2]))
        log.write("\n")


def read_input(file_name):
    try:
        file=open(file_name,"r")
    except IOError:
        error("There is a problema reading the {file_name} file".format(file_name=file_name))
        exit()


def main():
    
    workFlow=workflow.workFlow()
    
    write_log("log.file", workFlow)

    write_mdp(workFlow)
    write_run(run_file='runmd.sh', workdir='path/to/work/dir', program='path/to/gromacs', init_struct='initial_structure', topo='topology', mdp='path/to/mdp/directory', workflow=workFlow)
    write_submission(file_name="runjob", job_name="TEST1", job_file="runmd.sh")


def error(message):
    print("An error occurred.")
    print(message)
    print("See './__main__.py -h' for more instructions.")
    exit()
