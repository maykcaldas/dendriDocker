#!/bin/python3

import os
#import numpy as np
#import matplotlib as plot


ff=0
lf=1000

cases={
    "5-Fluorouracil" : [4, 5],
    "Carbamazepine" : [4],
    "Quercetin" : [0, 1, 2, 3],
    "Methotrexate" : [4],
    "SilybinA" : [2, 3, 4],
}
ligand={
    "5-Fluorouracil" : "2S04",
    "Carbamazepine" : "",
    "Quercetin" : "",
    "Methotrexate" : "",
    "SilybinA" : "",
}

#os.system("gmx514") GROMACS 5.1.4

print(cases)

################# make index ###############

# os.system("gmx make_ndx -f npt300_2.gro -o index.ndx << !\n\
# keep 0 \n\
# r CORE | r INTR | r TER \n\
# name 1 Dendrimer \n\
# r {ligand} \n\
# name 2 ligand \n\
# \n\
# q \n\
# !".format(ligand=ligand["5-Fluorouracil"]))

# print()


# ###################### RG ###################

# os.system("echo Dendrimer | gmx gyrate  -s md.tpr \
#                                         -f md.trr \
#                                         -n index.ndx \
#                                         -b {ff} \
#                                         -e {lf} \
#                                         -o gyrate.xvg".format(ff=ff, lf=lf))
# os.system("mv foo bar")

file = open("gyrate.xvg",'r')

Rg, Rg_x, Rg_y, Rg_z = 0,0,0,0
Rg2, Rg_x2, Rg_y2, Rg_z2 = 0,0,0,0
count = 0

for line in file:
    if line[0] not in ["@", "#"]:
        count+=1
        line_elements=line.strip().split()
        Rg+=float(line_elements[1])
        Rg2+=float(line_elements[1])*float(line_elements[1])

        Rg_x+=float(line_elements[2])
        Rg_x2+=float(line_elements[2])*float(line_elements[2])

        Rg_y+=float(line_elements[3])
        Rg_y2+=float(line_elements[3])*float(line_elements[3])

        Rg_z+=float(line_elements[4])
        Rg_z2+=float(line_elements[4])*float(line_elements[4])

Rgm   = Rg/count
Rg_dp = (Rg2/count - (Rg/count)**2 )**0.5
Rg_xm = Rg_x/count
Rg_xm_dp = (Rg_x2/count - (Rg_x/count)**2 )**0.5
Rg_ym = Rg_y/count
Rg_ym_dp = (Rg_y2/count - (Rg_y/count)**2 )**0.5
Rg_zm = Rg_z/count
Rg_zm_dp = (Rg_z2/count - (Rg_z/count)**2 )**0.5

file.close()

print("Rg:  ", Rgm, "+/-", Rg_dp)
print("Rgx: ", Rg_xm, "+/-", Rg_xm_dp)
print("Rgy: ", Rg_ym, "+/-", Rg_ym_dp)
print("Rgz: ", Rg_zm, "+/-", Rg_zm_dp)


# #calculate mean values

# ################### delta ###################

# # O ALGORITMO DE CALCULO DO MOMENTO DE INERCIA DO GROMACS TÁ BUGADO.
# # USAR MINHA ROTINA DE CÁLCULO DO DELTA

# os.system("echo Dendrimer | gmx gyrate  -moi \
#                                         -s md.tpr \
#                                         -f md.trr \
#                                         -n index.ndx \
#                                         -b {ff} \
#                                         -e {lf} \
#                                         -o moi.xvg".format(ff=ff, lf=lf))
# os.system("mv foo bar")

# #calculate mean values

# # ################ Create .gro ################

# # os.system("echo Dendrimer | $gmx514 trjconv -f md.xtc \
# #                                             -s md.tpr \
# #                                             -n index.ndx \
# #                                             -o ${dend}_G${i}_${j}.gro \
# #                                             -pbc mol")


# # os.system("mv ${dend}_G${i}_${j}.gro /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.")

# ################### RDF #####################

# os.system("gmx rdf  -f md.xtc \
#                     -s md.tpr \
#                     -n index.ndx \
#                     -b {ff} \
#                     -e {lf} \
#                     -ref Dendrimer \
#                     -sel ligand \
#                     -bin 0.002 \
#                     -selrpos res_com \
#                     -o rdf.xvg".format(ff=ff, lf=lf))

################ distances ##################

# #                       -cutoff 5 \

# os.system("gmx pairdist -f md.xtc \
#                         -s md.tpr \
#                         -n index.ndx \
#                         -b {ff} \
#                         -e {lf} \
#                         -selrpos res_com \
#                         -ref Dendrimer \
#                         -refgrouping mol \
#                         -sel ligand \
#                         -selgrouping res \
#                         -o distances.xvg".format(ff=ff, lf=lf))

file = open("distances.xvg",'r')

time=[]
ligands_n=[]
Rg=Rgm
for line in file:
    if line[0] not in ["@", "#"]:
        line_elements=line.strip().split()
        t = float(line_elements[0])
        count=0
        for k in line_elements:
            if float(k) < Rg:
                count+=1
        time.append(t)
        ligands_n.append(count)

file.close()

print(time)
print(ligands_n)

#plot ligands_n x time

################# clean #######################

os.system("rm \#*")