#!/bin/bash

# Here we define some important variables.

HERE=/scratch/60061a/maykcaldas/Dendrimer/Quercetin/G2/Acid
PROGRAM=/scratch/60061a/gromacs-2016.4-GITHUB/LINUX/bin/gmx_mpi
WORKDIR=${HERE}/tmp

TOPO=${WORKDIR}/topol.top

MDP=${HERE}/../../../../mdp

pwd

mkdir -p ${WORKDIR}
cp ${HERE}/topol.top ${WORKDIR}
cp ${HERE}/dock.gro ${WORKDIR}/dock.gro
cp ${HERE}/PAMAM_G2_Acid.itp ${WORKDIR}
cp ${HERE}/quercetin.itp ${WORKDIR}
cp ${HERE}/plumed.dat ${WORKDIR}
cd ${WORKDIR}


################### BOX ######################
#
#${PROGRAM} editconf \
#        -f ${HERE}/complex.gro \
#        -c \
#        -d 1.0 \
#        -bt cubic \
#        -o box1.gro
#
########## VACCUM MINIMIZATION ###############
#
#${PROGRAM} grompp \
#    -f ${MDP}/em.mdp \
#    -c ${WORKDIR}/box1.gro \
#    -p ${TOPO} \
#    -maxwarn 5 \
#    -o em1.tpr
#
#${PROGRAM} mdrun \
#    -s em1.tpr \
#    -deffnm em1
#
#
#################### BOX ######################
#
#${PROGRAM} editconf \
#	-f ${WORKDIR}/em1.gro \
#	-c \
#	-d 1.0 \
#	-bt cubic \
#	-o box2.gro
#
################### SOLVATE ####################
#
#${PROGRAM} solvate \
#	-cp ${WORKDIR}/box2.gro \
#	-cs spc216.gro \
#	-p ${TOPO} \
#	-o solv.gro
#	
##################### ION ######################
#
#${PROGRAM} grompp \
#	-f ${MDP}/ion.mdp \
#	-c ${WORKDIR}/solv.gro \
#	-p ${TOPO} \
#	-o ion.tpr
#
#
#echo SOL | ${PROGRAM} genion \
#	-s ${WORKDIR}/ion.tpr \
#	-p ${TOPO} \
#	-pname NA \
#	-nname CL \
#	-neutral \
#	-o ion.gro
#
########### ENERGY MINIMIZATION ##############

${PROGRAM} grompp \
    -f ${MDP}/em.mdp \
    -c ${WORKDIR}/dock.gro \
    -p ${TOPO} \
    -maxwarn 5 \
    -o em2.tpr

${PROGRAM} mdrun \
    -s em2.tpr \
    -deffnm em2

######### EQUILIBRATION 1: NVT50  ############

${PROGRAM} grompp \
    -f ${MDP}/nvt50.mdp \
    -c ${WORKDIR}/em2.gro \
    -p ${TOPO} \
    -maxwarn 2 \
    -o nvt50.tpr

${PROGRAM} mdrun \
    -s nvt50.tpr \
    -deffnm nvt50
    
######## EQUILIBRATION 1: NVT100  ############

${PROGRAM} grompp \
    -f ${MDP}/nvt100.mdp \
    -c ${WORKDIR}/nvt50.gro \
    -p ${TOPO} \
    -maxwarn 2 \
    -o nvt100.tpr

${PROGRAM} mdrun \
    -s nvt100.tpr \
    -deffnm nvt100

######### EQUILIBRATION 1: NVT150  ###########

${PROGRAM} grompp \
    -f ${MDP}/nvt150.mdp \
    -c ${WORKDIR}/nvt100.gro \
    -p ${TOPO} \
    -maxwarn 2 \
    -o nvt150.tpr

${PROGRAM} mdrun \
    -s nvt150.tpr \
    -deffnm nvt150

######### EQUILIBRATION 1: NVT200  ###########

${PROGRAM} grompp \
    -f ${MDP}/nvt200.mdp \
    -c ${WORKDIR}/nvt150.gro \
    -p ${TOPO} \
    -maxwarn 2 \
    -o nvt200.tpr

${PROGRAM} mdrun \
    -s nvt200.tpr \
    -deffnm nvt200

######### EQUILIBRATION 1: NVT250  ###########

${PROGRAM} grompp \
    -f ${MDP}/nvt250.mdp \
    -c ${WORKDIR}/nvt200.gro \
    -p ${TOPO} \
    -maxwarn 2 \
    -o nvt250.tpr

${PROGRAM} mdrun \
    -s nvt250.tpr \
    -deffnm nvt250
    
######### EQUILIBRATION 1: NVT300  ###########

${PROGRAM} grompp \
    -f ${MDP}/nvt300.mdp \
    -c ${WORKDIR}/nvt250.gro \
    -p ${TOPO} \
    -maxwarn 2 \
    -o nvt300.tpr

${PROGRAM} mdrun \
    -s nvt300.tpr \
    -deffnm nvt300
    
######### EQUILIBRATION 2: NPT300  ##############

${PROGRAM} grompp \
    -f ${MDP}/npt300.mdp \
    -c ${WORKDIR}/nvt300.gro \
    -p ${TOPO} \
    -maxwarn 5 \
    -o npt300.tpr

${PROGRAM} mdrun \
    -s npt300.tpr \
    -deffnm npt300

######### EQUILIBRATION 2: NPT400  ##############

${PROGRAM} grompp \
    -f ${MDP}/npt400.mdp \
    -c ${WORKDIR}/npt300.gro \
    -p ${TOPO} \
    -maxwarn 5 \
    -o npt400.tpr

${PROGRAM} mdrun \
    -s npt400.tpr \
    -deffnm npt400

######### EQUILIBRATION 2: NPT300  ##############

${PROGRAM} grompp \
    -f ${MDP}/npt300.mdp \
    -c ${WORKDIR}/npt400.gro \
    -p ${TOPO} \
    -maxwarn 5 \
    -o npt300_2.tpr

${PROGRAM} mdrun \
    -s npt300_2.tpr \
    -deffnm npt300_2

########## MOLECULAR DYNAMICS ###############

${PROGRAM} grompp \
    -f ${MDP}/md.mdp \
    -c ${WORKDIR}/npt300_2.gro \
    -p ${TOPO} \
    -maxwarn 5 \
    -o md.tpr

${PROGRAM} mdrun \
    -s md.tpr \
    -plumed plumed.dat \
    -cpi state.cpt \
    -cpo state.cpt \
    -deffnm md

#mkdir -p ${HERE}/em ${HERE}/nvt ${HERE}/npt ${HERE}/md ${HERE}/out
#
#mv em.* ${HERE}/em
#mv nvt.* ${HERE}/nvt
#mv npt.* ${HERE}/npt
#mv md.* ${HERE}/md
#mv *.* ${HERE}/out
