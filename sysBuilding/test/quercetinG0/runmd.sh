#!/bin/bash

# Here we define some important variables.

HERE=.
PROGRAM=/home/mayk/programs/gromacs-5.1.4/install/bin/gmx
WORKDIR=${HERE}

TOPO=topol.top
MDP=./mdp

cd ${HERE}


################### BOX ######################

${PROGRAM} editconf \
        -f $HERE/PAMAM_G0_Neutral.gro \
        -c \
        -d 1.0 \
        -bt cubic \
        -o box1.gro\

########## INSERT MOLECULES ##############
${PROGRAM} insert-molecules \
	 	 -f $WORKDIR/box1.gro \
	 	 -ci quercetin.pdb \
	 	 -nmol 8 \
	 	 -o box2\

################### BOX ######################

${PROGRAM} editconf \
        -f $HERE/box2.gro \
        -c \
        -d 0.2 \
        -bt cubic \
        -o box3.gro\

########## MINIMIZATION ###############

${PROGRAM} grompp \
    -f ${MDP}/em.mdp \
    -c $WORKDIR/box3.gro \
    -p ${TOPO} \
    -maxwarn 5 \
    -o em1.tpr\

${PROGRAM} mdrun \
    -deffnm em1\


################### BOX ######################

${PROGRAM} editconf \
        -f $HERE/em1.gro \
        -c \
        -d 0.2 \
        -bt cubic \
        -o box4.gro\

################### SOLVATE ####################

${PROGRAM} solvate \
	-cp $WORKDIR/box4.gro \
	-cs spc216.gro \
	-p ${TOPO} \
	-o solv.gro\

##################### ION ######################

${PROGRAM} grompp \
	-f ${MDP}/ion.mdp \
	-c $WORKDIR/solv.gro \
	-p ${TOPO} \
	-o ion.tpr\


echo SOL | ${PROGRAM} genion \
	-s ${WORKDIR}/ion.tpr \
	-p ${TOPO} \
	-pname NA+ \
	-nname CL- \
	-np  0 \
	-nn  0 \
	-o ion.gro\

########## MINIMIZATION ###############

${PROGRAM} grompp \
    -f ${MDP}/em.mdp \
    -c $WORKDIR/ion.gro \
    -p ${TOPO} \
    -maxwarn 5 \
    -o em2.tpr\

${PROGRAM} mdrun \
    -deffnm em2\


########## MOLECULAR DYNAMICS ###############

${PROGRAM} grompp \
         -f $MDP/dock.mdp \
         -c $WORKDIR/em2.gro \
         -p ${TOPO} \
         -maxwarn 2 \
         -o dock.tpr\

mpirun -n 8 $PROGRAM mdrun \
         -cpo dock.cpt \
         -deffnm dock\
         -plumed plumed.dat\

