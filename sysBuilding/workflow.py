#!/usr/bin/python3
#!-*- coding: utf8 -*-

'''
Software general documentation

gromacsBuilding was made using python 3.5.2

'''

def workFlow():
    workflow = [('BOX', 
        {'file_name': None},
        {'run_file': 'test/runmd.sh', 'init_struct': 'PAMAM_G0_Neutral.gro', 'd': '1.0', 'output': 'box'}
    ),
    # (SOLVATE, 
    #     {'file_name': None},
    #     {'run_file':'runmd.sh', 'system': 'box.gro', 'output': 'solv'}
    # ),
    ('ION', 
        {'file_name': None},
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
        {'run_file': 'runmd.sh', 'mdp': 'path/to/mdp', 'system': 'npt100', 'output': 'md', 'mpi': True, 'mpithreads': '8'}
    ),
    ]

    return workflow