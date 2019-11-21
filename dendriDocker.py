#!/usr/bin/python3
#!-*- coding: utf8 -*-

'''
Software general documentation

DendriBuilder was made using python 3.5.2

'''

from printer import printer
import argparse
import time
import os
import sysBuilding # Future plans to automatize the whole process.


def main():
    startTime=time.time()
    script_dir = os.path.dirname(__file__)
    logo=open(script_dir+'/asciilogo.txt','r')
    print(logo.read()[-285:])
    print("\nIt's beeing developed by the MSSM/LabMMol laboratory at IQ/UFRJ.")
    print("More about us in our webpage: https://labmmol.iq.ufrj.br/mssm/ \n")
    print("You need to have the topology and the coordinates files as well as the gromacs and plummed installed\n")
    print("in way to use this software.\n")
    print("Use -h to access the help menu.\n")

    parser=argparse.ArgumentParser(description="There is no description")
    # parser.add_argument('--coreSize', help="Number of atoms in the core block of the dendrimer")
    # parser.add_argument('--interSize', help="Number of atoms in the intermediary block of the dendrimer")
    # parser.add_argument('--terSize', help="Number of atoms in the terminal block of the dendrimer")
    parser.add_argument('--host', help="Name of the host molecule", default="Dendrimer")
    parser.add_argument('--ligand', help="Name of the molecule to be docked", default="Ligand")
    parser.add_argument('--nligand', help="Number of ligand molecules to be docked", default="10")
    parser.add_argument('--atomInDend', help="Number of atoms in your model of dendrimer", default=None)
    parser.add_argument('--atomInLigand', help="Number of atoms in your model of ligand", default=None)
    parser.add_argument('--ligandCoord', help="Path to the ligand coordinates file", default="ligand.gro")
    parser.add_argument('--dendCoord', help="Path to the dendrimer coordinates file", default="dend.gro")
    # parser.add_argument('--dendGeneration', help="Dendrimer generation", default="3")
    parser.add_argument('--ligandTop', help="Path to the ligand topology file", default="ligand.itp")
    parser.add_argument('--dendTop', help="Path to the dendrimer topology file", default="dend.itp")
    parser.add_argument('--FFPath', help="Path to the forcefield directory", default=None)
    parser.add_argument('--gmxPath', help="Path to the Gromacs binary file", default='/usr/local/gromacs')
    parser.add_argument('--mdpPath', help="Path to the mdp files directory", default='./mdp')
    parser.add_argument('--nameOut', help="Output name of the system", default='System')
    parser.add_argument('--topolOut', help="Output name of the topology file", default='topol.top')
    parser.add_argument('--dockOut', help="Output name of the docking file", default='plumed_dock.dat')
    parser.add_argument('--runOut', help="Output name of the run file", default='run.sh')
    parser.add_argument('--force', help="Force constant of the harmonic potential (or the slope in linear potential) in steered dynamics", default="500.0")
    parser.add_argument('--method', help="Docking method: harmonic, linear, shell, harmonicWall or linearWall", default="harmonic")
    
    args=vars(parser.parse_args())

    #Getting the atomInDend
    if args['dendCoord'][-4:] == '.gro':
        f=open(args['dendCoord'])
        f.readline()
        args['atomInDend']=f.readline().strip()
        f.close()

    noneInArgs=False
    for k in args:
        if args[k] is None:
            print("The {0:^15s} should be defined.".format(k))
            noneInArgs=True
    if noneInArgs == True:
        error("There were some undefined variables.")
    
    if not args['method'] in ['harmonic','linear','shell', 'harmonicWall', 'linearWall']:
        error("This is not a valid method.")

    print("The used list of parameters are:\n")
    for i in args:
        print("{0:15s}: {1}".format(i, args[i]))
    
    print("\n")    

    prt=printer()
    prt.printPlumedDock(args)
    prt.printTop(args)
    prt.printSetup(args)

    # I can create a workflow here and call the sysBuilding.
    # But I need to fix the sysBuilding directory structure
    # and include the sysBuilding inputs here in the dendriDocker.

    os.system('chmod u+x {0}'.format(args['runOut']))
    os.system('./{0}'.format(args['runOut']))
    os.system('rm -rf \#*')
    
    print("The run took: {0:.2g} seconds.\n".format(time.time()-startTime))


def error(message):
    print("An error occurred.")
    print(message)
    print("See './__main__.py -h' for more instructions.")
    exit()
