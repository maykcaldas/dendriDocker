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


def main():
    script_dir = os.path.dirname(__file__)
    logo=open(script_dir+'/asciilogo.txt','r')
    print(logo.read()[-285:])
    print("\nIt's beeing develop by the MSSM/LabMMol laboratory at IQ/UFRJ.")
    print("More about us in our webpage: https://labmmol.iq.ufrj.br/mssm/ \n")
    print("You need to have the topology and the coordinates files in way to use this software.")
    print("Use -h to access the help menu.\n")

    parser=argparse.ArgumentParser(description="There is no description")
    # parser.add_argument('--coreSize', help="Number of atoms in the core block of the dendrimer")
    # parser.add_argument('--interSize', help="Number of atoms in the intermediary block of the dendrimer")
    # parser.add_argument('--terSize', help="Number of atoms in the terminal block of the dendrimer")
    parser.add_argument('--ligand', help="PDB of the molecule to be docked", default="ligand.pdb")
    parser.add_argument('--nligand', help="Number of ligand molecules to be docked", default="10")
    parser.add_argument('--atomInDend', help="Number of atoms in your model dendrimer", default=None)
    parser.add_argument('--atomInLigand', help="Number of atoms in your model ligand", default=None)
    parser.add_argument('--ligandCoord', help="Path to the ligand coordinates file", default="ligand.gro")
    parser.add_argument('--dendCoord', help="Path to the dendrimer coordinates file", default="dend.gro")
    parser.add_argument('--ligandTop', help="Path to the ligand topology file", default="ligand.itp")
    parser.add_argument('--dendTop', help="Path to the dendrimer topology file", default="dend.itp")
    parser.add_argument('--FFPath', help="Path to the forcefield directory", default=None)
    parser.add_argument('--gmxPath', help="Path to the Gromacs binary file", default='/usr/local/gromacs')
    parser.add_argument('--mdpPath', help="Path to the mdp files directory", default='./mdp')
    parser.add_argument('--nameOut', help="Output name of the system", default='System')
    parser.add_argument('--topolOut', help="Output name of the topology file", default='topol.top')
    parser.add_argument('--dockOut', help="Output name of the docking file", default='plumed_dock.top')
    parser.add_argument('--runOut', help="Output name of the run file", default='run.sh')
    #parser.add_argument('--walls', help="Number of atoms in your model ligand")

    args=vars(parser.parse_args())

    noneInArgs=False
    for k in args:
        if args[k] is None:
            print("The {0} should be defined.".format(k))
            noneInArgs=True
    if noneInArgs == True:
        error("There were some undefined variables.")

    #coreSize=args["coreSize"]
    #interSize=args["interSize"]
    #terSize=args["terSize"]
    #ligand=args["ligand"]
    #nligand=args["nligand"]

    #os.system('ls -la')

    #print(args)
    prt=printer()
    prt.printPlumedDock(args)
    prt.printTop(args)
    prt.printSetup(args)


def error(message):
    print("An error occurred.")
    print(message)
    print("See './__main__.py -h' for more instructions.")
    exit()
