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
import sys

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

    workflow=[]

    for line in file:
        if line[:-1].strip() == "":
            continue
        elif line[:-1].strip()[0] == "#":
            continue
        elif line[:-1].strip() in ["BOX","INSERT","SOLV","ION","EM","NVT","NPT","MD"]:
            action=line[:-1].strip()
            
            op_dict={}
            input=file.readline()[:-1] #read first input options
            if input == "mdp_input":
                options=file.readline()[:-1].split(',')
                for option in options:
                    op_dict[option.split(':')[0].strip()]=option.split(':')[1].strip()
                mdp_options=op_dict
            elif input=="run_input":
                options=file.readline()[:-1].split(',')
                for option in options:
                    op_dict[option.split(':')[0].strip()]=option.split(':')[1].strip()
                run_options=op_dict
            else:
                error("It is expected to pass only the mdp_input and the run_input")

            op_dict={}
            input=file.readline()[:-1] #read second input options
            if input == "mdp_input":
                options=file.readline()[:-1].split(',')
                for option in options:
                    op_dict[option.split(':')[0].strip()]=option.split(':')[1].strip()
                mdp_options=op_dict
            elif input=="run_input":
                options=file.readline()[:-1].split(',')
                for option in options:
                    op_dict[option.split(':')[0].strip()]=option.split(':')[1].strip()
                run_options=op_dict
            else:
                error("It is expected to pass only the mdp_input and the run_input")
        else:
            error("There is an error in the input file")
        
        workflow.append((action,mdp_options,run_options))
        
    # for k in workflow:
    #     print(k)
    return workflow


def main():
    
    workFlow=read_input(sys.argv[1])

    # workFlow=create_workFlow()
    
    write_log("log.file", workFlow)

    write_mdp(workFlow)
    write_run(run_file='runmd.sh', workdir='path/to/work/dir', program='path/to/gromacs', init_struct='initial_structure', topo='topology', mdp='path/to/mdp/directory', workflow=workFlow)
    write_submission(file_name="runjob", job_name="TEST1", job_file="runmd.sh")

def create_workFlow():
    workflow = [('BOX', 
        {'file_name': None},
        {'run_file': 'test/runmd.sh', 'init_struct': 'PAMAM_G0_Neutral.gro', 'd': '1.0', 'output': 'box'}
    ),
    # (SOLVATE, 
    #     {'file_name': None},
    #     {'run_file':'runmd.sh', 'system': 'box.gro', 'output': 'solv'}
    # ),
    ('ION', 
        {'file_name': 'None'},
        {'run_file': 'runmd.sh', 'system': 'box', 'output': 'ion', 'neutral':True, 'na':'0', 'cl':'0'}
    ),
    ('EM', 
        {'file_name': 'em.mdp', 'emtol': '100.0', 'nsteps': '5000', 'emstep': '0.001'}, 
        {'run_file': 'runmd.sh', 'system': 'ion', 'output': 'em1'}
    ),
    ('NVT', 
        {'file_name': 'nvt100.mdp', 'nsteps': '25000000', 'dt': '0.002', 'cont':'no', 'temp':'100'}, 
        {'run_file': 'runmd.sh', 'mdp': 'path/to/mdp', 'system': 'em1', 'output': 'nvt100', 'mpi': True, 'mpithreads': '8'}
    ),
    ('NPT', 
        {'file_name': 'npt100.mdp', 'nsteps': '100000', 'dt': '0.002', 'cont':'yes', 'temp':'100', 'press': '1'}, 
        {'run_file': 'runmd.sh', 'mdp': 'path/to/mdp', 'system': 'nvt100', 'output': 'npt100', 'mpi': True, 'mpithreads': '8'}
    ),
    ('MD', 
        {'file_name': 'md.mdp', 'nsteps': '100000', 'dt': '0.002', 'cont':'yes', 'temp':'298', 'press': '1'}, 
        {'run_file': 'runmd.sh', 'mdp': 'path/to/mdp', 'system': 'npt100', 'output': 'md', 'mpi': True, 'mpithreads': '8', 'plumed': True, 'plumed_file': 'plumed.dat'}
    ),
    ]

    return workflow

def error(message):
    print("An error occurred.")
    print(message)
    print("See './__main__.py -h' for more instructions.")
    exit()
