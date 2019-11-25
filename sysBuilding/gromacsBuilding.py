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
