#!/usr/bin/python3
#!-*- coding: utf8 -*-

'''
Software general documentation

gromacsBuilding was made using python 3.5.2

'''

import os

def write_submission(file_name, job_name, run_file):
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
    job.write('{run_file}'.format(run_file=run_file))


def write_run_box(run_file, mdp_file, init_struct, d, output):
    run_file.write('################### BOX ######################\n')
    run_file.write('\n')
    run_file.write('${PROGRAM} editconf \\\n')
    run_file.write('        -f ${HERE}/{init_struct} \\\n'.format(HERE="HERE",init_struct=init_struct))
    run_file.write('        -c \\\n')
#    run_file.write('        -d {d} \\\n'.format(d=d))
    run_file.write('        -bt cubic \\\n')
    run_file.write('        -o {output}.gro\\\n'.format(output=output))
    run_file.write('\n')


def write_run_insert_molecules(run_file, mdp_file, system, insert, nmol, output):
    run_file.write('########## INSERT MOLECULES ##############\n')
    run_file.write('${PROGRAM} insert-molecules \\\n')
    run_file.write('\t \t -f ${WORKDIR}/{system} \\\n'.format(WORKDIR="WORKDIR", system=system))
    run_file.write('\t \t -ci {insert} \\\n'.format(insert=insert))
    run_file.write('\t \t -nmol {nmol} \\\n'.format(nmol=nmol))
    run_file.write('\t \t -o {output}\\\n'.format(output=output))
    run_file.write('\n')


def write_run_solvate(run_file, mdp_file, system, output):
    run_file.write('################### SOLVATE ####################\n')
    run_file.write('\n')
    run_file.write('${PROGRAM} solvate \\\n')
    run_file.write('	-cp ${WORKDIR}/{system} \\\n'.format(WORKDIR="WORKDIR", system=system))
    run_file.write('	-cs spc216.gro \\\n')
    run_file.write('	-p ${TOPO} \\\n')
    run_file.write('	-o {output}.gro\\\n'.format(output=output))
    run_file.write('\n')


def write_run_ion(run_file, mdp_file, system, output, neutral='True', na=0, cl=0):
    run_file.write('##################### ION ######################\n')
    run_file.write('\n')
    run_file.write('${PROGRAM} grompp \\\n')
    run_file.write('    -f ${MDP}/{mdp} \\\n'.format(MDP="MDP", mdp=mdp_file))
    run_file.write('	-c ${WORKDIR}/{system} \\\n'.format(WORKDIR="WORKDIR", system=system))
    run_file.write('	-p ${TOPO} \\\n')
    run_file.write('	-maxwarn 2 \\\n')
    run_file.write('	-o {output}.tpr\\\n'.format(output=output))
    run_file.write('\n')
    run_file.write('\n')
    
    run_file.write('echo SOL | ${PROGRAM} genion \\\n')
    run_file.write('	-s ${WORKDIR}/ion.tpr \\\n')
    run_file.write('	-p ${TOPO} \\\n')
    run_file.write('	-pname NA+ \\\n')
    run_file.write('	-nname CL- \\\n')
    if neutral == 'True':
        run_file.write('	-neutral \\\n')
    else:
        run_file.write('	-np  {na} \\\n'.format(na=na))
        run_file.write('	-nn  {cl} \\\n'.format(cl=cl))
    run_file.write('	-o {output}.gro\\\n'.format(output=output))
    run_file.write('\n')


def write_run_em(run_file, mdp_file, system, output):
    run_file.write('########## MINIMIZATION ###############\n')
    run_file.write('\n')
    run_file.write('${PROGRAM} grompp \\\n')
    run_file.write('    -f ${MDP}/{mdp} \\\n'.format(MDP="MDP", mdp=mdp_file))
    run_file.write('    -c ${WORKDIR}/{system} \\\n'.format(WORKDIR="WORKDIR",system=system))
    run_file.write('    -p ${TOPO} \\\n')
    run_file.write('    -maxwarn 5 \\\n')
    run_file.write('    -o {output}.tpr\\\n'.format(output=output))
    run_file.write('\n')
    
    run_file.write('${PROGRAM} mdrun \\\n')
    run_file.write('    -deffnm {output}\\\n'.format(output=output))
    run_file.write('\n')
    run_file.write('\n')


def write_run_nvt(run_file, mdp_file, system, output, mpi='True', mpithreads=8):
    run_file.write('######### EQUILIBRATION: NVT  ############\n')
    run_file.write('\n')
    run_file.write('${PROGRAM} grompp \\\n')
    run_file.write('         -f ${MDP}/{mdp} \\\n'.format(MDP="MDP", mdp=mdp_file))
    run_file.write('         -c ${WORKDIR}/{system} \\\n'.format(WORKDIR="WORKDIR",system=system))
    run_file.write('         -p ${TOPO} \\\n')
    run_file.write('         -maxwarn 2 \\\n')
    run_file.write('         -o {output}.tpr\\\n'.format(output=output))
    run_file.write('\n')
    if mpi=='True':
        run_file.write('mpirun -n {mpithreads} ${PROGRAM} mdrun \\\n'.format(mpithreads=mpithreads, PROGRAM="PROGRAM"))
    else:
        run_file.write('${PROGRAM} mdrun \\\n')
    run_file.write('         -deffnm {output}\\\n'.format(output=output))
    run_file.write('\n')


def write_run_npt(run_file, mdp_file, system, output, mpi='True', mpithreads=8):
    run_file.write('######## EQUILIBRATION: NPT  ##############\n')
    run_file.write('\n')
    run_file.write('${PROGRAM} grompp \\\n')
    run_file.write('         -f ${MDP}/{mdp} \\\n'.format(MDP="MDP", mdp=mdp_file))
    run_file.write('         -c ${WORKDIR}/{system} \\\n'.format(WORKDIR="WORKDIR",system=system))
    run_file.write('         -p ${TOPO} \\\n')
    run_file.write('         -maxwarn 2 \\\n')
    run_file.write('         -o {output}.tpr\\\n'.format(output=output))
    run_file.write('\n')
    if mpi=='True':
        run_file.write('mpirun -n {mpithreads} ${PROGRAM} mdrun \\\n'.format(mpithreads=mpithreads, PROGRAM="PROGRAM"))
    else:
        run_file.write('${PROGRAM} mdrun \\\n')
    run_file.write('         -deffnm {output}\\\n'.format(output=output))
    run_file.write('\n')


def write_run_md(run_file, mdp_file, system, output, mpi='True', mpithreads=8, plumed='True', plumed_file='plumed.dat'):
    run_file.write('########## MOLECULAR DYNAMICS ###############\n')
    run_file.write('\n')
    run_file.write('${PROGRAM} grompp \\\n')
    run_file.write('         -f ${MDP}/{mdp} \\\n'.format(MDP="MDP", mdp=mdp_file))
    run_file.write('         -c ${WORKDIR}/{system} \\\n'.format(WORKDIR="WORKDIR",system=system))
    run_file.write('         -p ${TOPO} \\\n')
    run_file.write('         -maxwarn 2 \\\n')
    run_file.write('         -o {output}.tpr\\\n'.format(output=output))
    run_file.write('\n')
    if mpi=='True':
        run_file.write('mpirun -n {mpithreads} ${PROGRAM} mdrun \\\n'.format(mpithreads=mpithreads, PROGRAM="PROGRAM"))
    else:
        run_file.write('${PROGRAM} mdrun \\\n')
    # run_file.write('         -cpo {output}.cpt \\\n'.format(output=output))
    run_file.write('         -deffnm {output}\\\n'.format(output=output))
    if plumed=='True':
        run_file.write('         -plumed {plumed_file}\\\n'.format(plumed_file=plumed_file))
    run_file.write('\n')


def write_run(run_file, workdir, program, init_struct, topo, mdpPath, workflow):
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
    run_file.write('WORKDIR=${HERE}\n')
    run_file.write('\n')
    run_file.write('TOPO={topo}\n'.format(topo=topo))
    run_file.write('MDP={mdpPath}\n'.format(mdpPath=mdpPath))
    run_file.write('\n')
    run_file.write('cd ${WORKDIR}\n')
    run_file.write('\n')
    run_file.write('\n')
    
    for step in workflow:
        if step[0] == "BOX":
            box_mdp_file=step[2]["mdp_file"]
            box_init_structure=step[2]["init_struct"]
            box_d=step[2]["d"]
            box_output=step[2]["output"]

            write_run_box(run_file, box_mdp_file, box_init_structure, box_d, box_output)
        elif step[0] == "SOLV":
            solv_mdp_file=step[2]["mdp_file"]
            solv_system=step[2]["system"]
            solv_output=step[2]["output"]

            write_run_solvate(run_file, solv_mdp_file, solv_system, solv_output)
        elif step[0] == "EM":
            em_mdp_file=step[2]["mdp_file"]
            em_system=step[2]["system"]
            em_output=step[2]["output"]
            
            write_run_em(run_file, em_mdp_file, em_system, em_output)
        elif step[0] == "ION":
            ion_mdp_file=step[2]["mdp_file"]
            ion_system=step[2]["system"]
            ion_output=step[2]["output"]
            ion_neutral=step[2]["neutral"]
            ion_na=step[2]["na"]
            ion_cl=step[2]["cl"]

            write_run_ion(run_file, ion_mdp_file, ion_system, ion_output, ion_neutral, ion_na, ion_cl)
        elif step[0] == "INSERT":
            insert_mdp_file=step[2]["mdp_file"]
            insert_system=step[2]["system"]
            insert_insert=step[2]["insert"]
            insert_nmol=step[2]["nmol"]
            insert_output=step[2]["output"]

            write_run_insert_molecules(run_file, insert_mdp_file, insert_system, insert_insert, insert_nmol, insert_output)
        elif step[0] == "SD":
            pass
        elif step[0] == "NVT":
            nvt_mdp_file=step[2]["mdp_file"]
            nvt_system=step[2]["system"]
            nvt_output=step[2]["output"]
            nvt_mpi=step[2]["mpi"]
            nvt_mpiThreads=step[2]["mpithreads"]
            
            write_run_nvt(run_file, nvt_mdp_file, nvt_system, nvt_output, nvt_mpi, nvt_mpiThreads)
        elif step[0] == "NPT":
            npt_mdp_file=step[2]["mdp_file"]
            npt_system=step[2]["system"]
            npt_output=step[2]["output"]
            npt_mpi=step[2]["mpi"]
            npt_mpiThreads=step[2]["mpithreads"]
            
            write_run_npt(run_file, npt_mdp_file, npt_system, npt_output, npt_mpi, npt_mpiThreads)
        elif step[0] == "MD":
            md_mdp_file=step[2]["mdp_file"]
            md_system=step[2]["system"]
            md_output=step[2]["output"]
            md_mpi=step[2]["mpi"]
            md_mpiThreads=step[2]["mpithreads"]
            md_plumed=step[2]["plumed"]
            md_plumed_file=step[2]["plumed_file"]
            
            write_run_md(run_file, md_mdp_file, md_system, md_output, md_mpi, md_mpiThreads, md_plumed, md_plumed_file)
        else:
            error("There is not such run process implemented. Please, check your input or contact the developers.")


def error(message):
    print("An error occurred.")
    print(message)
    print("See './__main__.py -h' for more instructions.")
    exit()
