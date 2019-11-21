#!/bin/bash

# Here we define some important variables.

HERE=path/to/work/dir
PROGRAM=path/to/gromacs
WORKDIR=${HERE}/tmp

TOPO=topology
MDP=path/to/mdp/directory

cd ${HERE}

mkdir -p ${WORKDIR}
cp ${HERE}/topol.top ${WORKDIR}
cp ${HERE}/dock.gro ${WORKDIR}/dock.gro
cp ${HERE}/PAMAM_G0_Neutral.itp ${WORKDIR}
cp ${HERE}/quercetin.itp ${WORKDIR}
#cp ${HERE}/plumed.dat ${WORKDIR}
cd ${WORKDIR}


################### BOX ######################

${PROGRAM} editconf 
        -f $HERE/PAMAM_G0_Neutral.gro.gro 
        -c 
        -d 1.0 
        -bt cubic 
        -o box.gro

##################### ION ######################

${PROGRAM} grompp 
	-f ${MDP}/ion.mdp 
	-c $WORKDIR/box.gro 
	-p ${TOPO} 
	-o ion.tpr


echo SOL | ${PROGRAM} genion 
	-s ${WORKDIR}/ion.tpr 
	-p ${TOPO} 
	-pname NA 
	-nname CL 
	-neutral 
	-o ion.gro

########## MINIMIZATION ###############

${PROGRAM} grompp 
    -f ${MDP}/em.mdp 
    -c $WORKDIR/ion.gro 
    -p ${TOPO} 
    -maxwarn 5 
    -o em1.tpr

${PROGRAM} mdrun_file 
    -deffnm em1


######### EQUILIBRATION: NVT  ############

${PROGRAM} grompp 
         -f $MDP/path/to/mdp 
         -c $WORKDIR/em1.gro 
         -p ${TOPO} 
         -maxwarn 2 
         -o nvt100.tpr

mpirun -n 8 $PROGRAM mdrun 
         -deffnm nvt100

######## EQUILIBRATION: NPT  ##############

${PROGRAM} grompp 
         -f $MDP/path/to/mdp.mdp 
         -c $WORKDIR/nvt100.gro 
         -p ${TOPO} 
         -maxwarn 2 
         -o npt100.tpr

mpirun -n 8 $PROGRAM mdrun 
         -deffnm npt100

########## MOLECULAR DYNAMICS ###############

${PROGRAM} grompp 
         -f $MDP/path/to/mdp.mdp 
         -c $WORKDIR/npt100.gro 
         -p ${TOPO} 
         -maxwarn 2 
         -o md.tpr

mpirun -n 8 $PROGRAM mdrun 
         -cpo md.cpt 
         -deffnm md

