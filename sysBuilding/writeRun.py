#!/usr/bin/python3
#!-*- coding: utf8 -*-

import os

'''
Software general documentation

gromacsBuilding was made using python 3.5.2

'''

def write_run_box(run_file, init_struct, d, output):
    run_file.write('################### BOX ######################\n')
    run_file.write('\n')
    run_file.write('${PROGRAM} editconf \n')
    run_file.write('        -f ${HERE}/{init_struct}.gro \n'.format(HERE="HERE",init_struct=init_struct))
    run_file.write('        -c \n')
    run_file.write('        -d {d} \n'.format(d=d))
    run_file.write('        -bt cubic \n')
    run_file.write('        -o {output}.gro\n'.format(output=output))
    run_file.write('\n')


def write_insert_molecules(run_file, system, insert, nmol, output):
    run_file.write('### INSERT MOLECULES ###\n')
    run_file.write('${gmx} insert-molecules \\\n')
    run_file.write('\t \t -f {system} \\\n'.format(system=system))
    run_file.write('\t \t -ci {} \\\n'.format(insert=insert))
    run_file.write('\t \t -nmol {} \\\n'.format(nmol=nmol))
    run_file.write('\t \t -o box.gro\n')
    run_file.write('\n')


def write_run_solvate(run_file, system, output):
    run_file.write('################### SOLVATE ####################\n')
    run_file.write('\n')
    run_file.write('$PROGRAM} solvate \\\n')
    run_file.write('	-cp ${WORKDIR}/{system}.gro \\\n'.format(WORKDIR="WORKDIR",system=system))
    run_file.write('	-cs spc216.gro \\\n')
    run_file.write('	-p ${TOPO} \\\n')
    run_file.write('	-o {output}.gro\n'.format(output=output))
    run_file.write('\n')


def write_run_ion(run_file, system, output, neutral=True, na=0, cl=0):
    run_file.write('##################### ION ######################\n')
    run_file.write('\n')
    run_file.write('${PROGRAM} grompp \n')
    run_file.write('	-f ${MDP}/ion.mdp \n')
    run_file.write('	-c ${WORKDIR}/{system}.gro \n'.format(WORKDIR="WORKDIR", system=system))
    run_file.write('	-p ${TOPO} \n')
    run_file.write('	-o {output}.tpr\n'.format(output=output))
    run_file.write('\n')
    run_file.write('\n')
    
    run_file.write('echo SOL | ${PROGRAM} genion \n')
    run_file.write('	-s ${WORKDIR}/ion.tpr \n')
    run_file.write('	-p ${TOPO} \n')
    run_file.write('	-pname NA \n')
    run_file.write('	-nname CL \n')
    if neutral == True:
        run_file.write('	-neutral \n')
    else:
        run_file.write('	-np  {na} \n'.format(na=na))
        run_file.write('	-nn  {cl} \n'.format(cl=cl))
    run_file.write('	-o {output}.gro\n'.format(output=output))
    run_file.write('\n')


def write_run_em(run_file, system, output):
    run_file.write('########## MINIMIZATION ###############\n')
    run_file.write('\n')
    run_file.write('${PROGRAM} grompp \n')
    run_file.write('    -f ${MDP}/em.mdp \n')
    run_file.write('    -c ${WORKDIR}/{system}.gro \n'.format(WORKDIR="WORKDIR",system=system))
    run_file.write('    -p ${TOPO} \n')
    run_file.write('    -maxwarn 5 \n')
    run_file.write('    -o {output}.tpr\n'.format(output=output))
    run_file.write('\n')
    
    run_file.write('${PROGRAM} mdrun_file \n')
    run_file.write('    -deffnm {output}\n'.format(output=output))
    run_file.write('\n')
    run_file.write('\n')


def write_run_nvt(run_file, mdp, system, output, mpi=True, mpithreads=8):
    run_file.write('######### EQUILIBRATION: NVT  ############\n')
    run_file.write('\n')
    run_file.write('${PROGRAM} grompp \n')
    run_file.write('         -f ${MDP}/{mdp} \n'.format(MDP="MDP",mdp=mdp))
    run_file.write('         -c ${WORKDIR}/{system}.gro \n'.format(WORKDIR="WORKDIR",system=system))
    run_file.write('         -p ${TOPO} \n')
    run_file.write('         -maxwarn 2 \n')
    run_file.write('         -o {output}.tpr\n'.format(output=output))
    run_file.write('\n')
    if mpi==True:
        run_file.write('mpirun -n {mpithreads} ${PROGRAM} mdrun \n'.format(mpithreads=mpithreads, PROGRAM="PROGRAM"))
    else:
        run_file.write('${PROGRAM} mdrun \n')
    run_file.write('         -deffnm {output}\n'.format(output=output))
    run_file.write('\n')


def write_run_npt(run_file, mdp, system, output, mpi=True, mpithreads=8):
    run_file.write('######## EQUILIBRATION: NPT  ##############\n')
    run_file.write('\n')
    run_file.write('${PROGRAM} grompp \n')
    run_file.write('         -f ${MDP}/{mdp}.mdp \n'.format(MDP="MDP",mdp=mdp))
    run_file.write('         -c ${WORKDIR}/{system}.gro \n'.format(WORKDIR="WORKDIR",system=system))
    run_file.write('         -p ${TOPO} \n')
    run_file.write('         -maxwarn 2 \n')
    run_file.write('         -o {output}.tpr\n'.format(output=output))
    run_file.write('\n')
    if mpi==True:
        run_file.write('mpirun -n {mpithreads} ${PROGRAM} mdrun \n'.format(mpithreads=mpithreads, PROGRAM="PROGRAM"))
    else:
        run_file.write('${PROGRAM} mdrun \n')
    # run_file.write('         -s npt.tpr \n')
    run_file.write('         -deffnm {output}\n'.format(output=output))
    run_file.write('\n')


def write_run_md(run_file, mdp, system, output, mpi=True, mpithreads=8, plumed=True, plumed_file='plumed.dat'):
    run_file.write('########## MOLECULAR DYNAMICS ###############\n')
    run_file.write('\n')
    run_file.write('${PROGRAM} grompp \n')
    run_file.write('         -f ${MDP}/{mdp}.mdp \n'.format(MDP="MDP",mdp=mdp))
    run_file.write('         -c ${WORKDIR}/{system}.gro \n'.format(WORKDIR="WORKDIR",system=system))
    run_file.write('         -p ${TOPO} \n')
    run_file.write('         -maxwarn 2 \n')
    run_file.write('         -o {output}.tpr\n'.format(output=output))
    run_file.write('\n')
    if mpi==True:
        run_file.write('mpirun -n {mpithreads} ${PROGRAM} mdrun \n'.format(mpithreads=mpithreads, PROGRAM="PROGRAM"))
    else:
        run_file.write('${PROGRAM} mdrun \n')
    # run_file.write('        -s md.tpr \n')
    # run_file.write('         -cpi state.cpt \n')
    run_file.write('         -cpo {output}.cpt \n'.format(output=output))
    run_file.write('         -deffnm {output}\n'.format(output=output))
    if plumed=True:
        run_file.write('         -plumed {plumed_file}\n'.format(plumed_file=plumed_file))
    run_file.write('\n')


def write_run(run_file, workdir, program, init_struct, topo, mdp, workflow):
    if os.path.isfile(run_file):
        print('!!!Backing up the existing run file!!!')
        os.rename(run_file, 'bck.'+run_file)
    run_file=open(run_file,'a')
    
    run_file.write('#!/bin/bash\n')
    run_file.write('\n')
    run_file.write('# Here we define some important variables.\n')
    run_file.write('\n')
    run_file.write('HERE={workdir}\n'.format(workdir=workdir))
    run_file.write('PROGRAM={program}\n'.format(program=program))
    run_file.write('WORKDIR=${HERE}/tmp\n')
    run_file.write('\n')
    run_file.write('TOPO={topo}\n'.format(topo=topo))
    run_file.write('MDP={mdp}\n'.format(mdp=mdp))
    run_file.write('\n')
    run_file.write('cd ${HERE}\n')
    run_file.write('\n')
    run_file.write('mkdir -p ${WORKDIR}\n')
    run_file.write('cp ${HERE}/topol.top ${WORKDIR}\n')
    run_file.write('cp ${HERE}/dock.gro ${WORKDIR}/dock.gro\n')
    run_file.write('cp ${HERE}/PAMAM_G0_Neutral.itp ${WORKDIR}\n')
    run_file.write('cp ${HERE}/quercetin.itp ${WORKDIR}\n')
    run_file.write('#cp ${HERE}/plumed.dat ${WORKDIR}\n')
    run_file.write('cd ${WORKDIR}\n')
    run_file.write('\n')
    run_file.write('\n')
    
    for step in workflow:
        if step[0] == "BOX":
            box_run_file=step[2]["run_file"]
            box_init_structure=step[2]["init_struct"]
            box_d=step[2]["d"]
            box_output=step[2]["output"]

            write_run_box(run_file, box_init_structure, box_d, box_output)
        elif step[0] == "SOLV":
            pass
        elif step[0] == "EM":
            em_run_file=step[2]["run_file"]
            em_system=step[2]["system"]
            em_output=step[2]["output"]
            
            write_run_em(run_file, em_system, em_output)
        elif step[0] == "ION":
            ion_run_file=step[2]["run_file"]
            ion_system=step[2]["system"]
            ion_output=step[2]["output"]
            ion_neutral=step[2]["neutral"]
            ion_na=step[2]["na"]
            ion_cl=step[2]["cl"]

            write_run_ion(run_file, ion_system, ion_output, ion_neutral, ion_na, ion_cl)
        elif step[0] == "SD":
            pass
        elif step[0] == "NVT":
            nvt_run_file=step[2]["run_file"]
            nvt_mdp=step[2]["mdp"]
            nvt_system=step[2]["system"]
            nvt_output=step[2]["output"]
            nvt_mpi=step[2]["mpi"]
            nvt_mpiThreads=step[2]["mpithreads"]
            
            write_run_nvt(run_file, nvt_mdp, nvt_system, nvt_output, nvt_mpi, nvt_mpiThreads)
        elif step[0] == "NPT":
            npt_run_file=step[2]["run_file"]
            npt_mdp=step[2]["mdp"]
            npt_system=step[2]["system"]
            npt_output=step[2]["output"]
            npt_mpi=step[2]["mpi"]
            npt_mpiThreads=step[2]["mpithreads"]
            
            write_run_npt(run_file, npt_mdp, npt_system, npt_output, npt_mpi, npt_mpiThreads)
        
        elif step[0] == "MD":
            md_run_file=step[2]["run_file"]
            md_mdp=step[2]["mdp"]
            md_system=step[2]["system"]
            md_output=step[2]["output"]
            md_mpi=step[2]["mpi"]
            md_mpiThreads=step[2]["mpithreads"]
            
            write_run_md(run_file, md_mdp, md_system, md_output, md_mpi, md_mpiThreads)
        
        else:
            error("There is not such run process implemented. Please, check your input or contact the developers.")
