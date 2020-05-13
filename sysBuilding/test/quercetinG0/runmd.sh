#!/bin/bash

# Here we define some important variables.

HERE=.
PROGRAM=/home/kyam/Programs/gromacs-2019.4/install/bin/gmx_mpi
WORKDIR=${HERE}

TOPO=topol.top
MDP=./mdp

cd ${HERE}


########## MINIMIZATION ###############

${PROGRAM} grompp \
    -f $MDP/em.mdp \
    -c $WORKDIR/dock.gro \
    -p ${TOPO} \
    -maxwarn 5 \
    -o em.tpr\

${PROGRAM} mdrun \
    -deffnm em\


######### EQUILIBRATION: NVT  ############

${PROGRAM} grompp \
         -f $MDP/nvt100.mdp \
         -c $WORKDIR/em \
         -p ${TOPO} \
         -maxwarn 2 \
         -o nvt100.tpr\

${PROGRAM} mdrun \
         -v \
         -deffnm nvt100\

######### EQUILIBRATION: NVT  ############

${PROGRAM} grompp \
         -f $MDP/nvt200.mdp \
         -c $WORKDIR/nvt100.gro \
         -p ${TOPO} \
         -maxwarn 2 \
         -o nvt200.tpr\

${PROGRAM} mdrun \
         -v \
         -deffnm nvt200\

######### EQUILIBRATION: NVT  ############

${PROGRAM} grompp \
         -f $MDP/nvt300.mdp \
         -c $WORKDIR/nvt200.gro \
         -p ${TOPO} \
         -maxwarn 2 \
         -o nvt300.tpr\

${PROGRAM} mdrun \
         -v \
         -deffnm nvt300\

######## EQUILIBRATION: NPT  ##############

${PROGRAM} grompp \
         -f $MDP/npt300.mdp \
         -c $WORKDIR/nvt300 \
         -p ${TOPO} \
         -maxwarn 2 \
         -o npt300.tpr\

${PROGRAM} mdrun \
         -v \
         -deffnm npt300\

########## MOLECULAR DYNAMICS ###############

${PROGRAM} grompp \
         -f $MDP/md.mdp \
         -c $WORKDIR/npt300.gro \
         -p ${TOPO} \
         -maxwarn 2 \
         -o md.tpr\

mpirun -n 2 $PROGRAM mdrun \
         -v \
         -deffnm md\

