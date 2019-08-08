#!/bin/python3

import os
import numpy as np
import matplotlib.pyplot as plt

#####   creating .gro  ###############


# # os.system("mv ${dend}_G${i}_${j}.gro /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.")


def index(ligand_case):
    os.system("gmx make_ndx -f npt300_2.gro -o index.ndx << !\n\
    keep 0 \n\
    r CORE | r INTR | r TER \n\
    name 1 Dendrimer \n\
    r {ligand} \n\
    name 2 ligand \n\
    \n\
    q \n\
    !".format(ligand=ligand_case))


def calculate_gyrate(outputName="gyrate.xvg", ff=0, lf=0):
    os.system("echo Dendrimer | gmx gyrate  -s md.tpr \
                                            -f md.xtc \
                                            -n index.ndx \
                                            -b {ff} \
                                            -e {lf} \
                                            -o {output}".format(output=outputName, ff=ff, lf=lf))
    file = open(outputName,'r')

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
    Rgm_dp = (Rg2/count - (Rg/count)**2 )**0.5
    Rg_xm = Rg_x/count
    Rg_xm_dp = (Rg_x2/count - (Rg_x/count)**2 )**0.5
    Rg_ym = Rg_y/count
    Rg_ym_dp = (Rg_y2/count - (Rg_y/count)**2 )**0.5
    Rg_zm = Rg_z/count
    Rg_zm_dp = (Rg_z2/count - (Rg_z/count)**2 )**0.5

    file.close()

    return ((Rgm, Rgm_dp),(Rg_xm, Rg_xm_dp),(Rg_ym, Rg_ym_dp),(Rg_zm, Rg_zm_dp))

def calculate_delta(outputName="output.gro", ff=0, lf=0):
    # O ALGORITMO DE CALCULO DO MOMENTO DE INERCIA DO GROMACS TÁ BUGADO.
    # USAR MINHA ROTINA DE CÁLCULO DO DELTA
    # Ela precisa ler um .gro. Não sei se eu quero fazer isso

    os.system("echo Dendrimer | gmx trjconv -f md.xtc \
                                            -s md.tpr \
                                            -b {ff} \
                                            -e {lf} \
                                            -n index.ndx \
                                            -pbc nojump\
                                            -o {output}".format(output=outputName, ff=ff, lf=lf))

    file=open(outputName, 'r')
    gro=file.readlines()
    file.close()


    

def rdf(outputName="rdf.xvg", ff=0, lf=0):
    os.system("gmx rdf  -f md.xtc \
                        -s md.tpr \
                        -n index.ndx \
                        -b {ff} \
                        -e {lf} \
                        -ref Dendrimer \
                        -sel ligand \
                        -bin 0.002 \
                        -selrpos res_com \
                        -o {output}".format(output=outputName, ff=ff, lf=lf))

    # Plot rdf

def distance(outputName="distances.xvg", ff=0, lf=0, Rgm=0):
    os.system("gmx pairdist -f md.xtc \
                            -s md.tpr \
                            -n index.ndx \
                            -b {ff} \
                            -e {lf} \
                            -selrpos res_com \
                            -ref Dendrimer \
                            -refgrouping mol \
                            -sel ligand \
                            -selgrouping res \
                            -o {output}".format(output=outputName, ff=ff, lf=lf))

    file = open(outputName,'r')
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

    return time, ligands_n


def clean(list):
    trash=""
    for c in list:
        trash+=c+" "
    os.system("rm {}".format(trash))

def main():
    ff=0
    lf=50000
    root="/home/mayk/Documents/Labmmol/Dendrimer/dendriDocker/validation"

    cases={
        "5-Fluorouracil" : [4, 5],
        "Carbamazepine" : [4],
        "Quercetin" : [0, 1],
        # "Methotrexate" : [4],
        # "SilybinA" : [2, 3, 4],
    }
    ligand={
        "5-Fluorouracil" : "2S04",
        "Carbamazepine" : "BDRM",
        "Quercetin" : "VV98",
        "Methotrexate" : "6QRE",
        "SilybinA" : "SYLI",
    }

    #os.system("gmx514") #Sourcing GROMACS 5.1.4

    for case in cases:
        for G in cases[case]:
            if case == "Carbamazepine":
                systems=["Acid/F_100", "Acid/F_500", "Neutral"]
            elif case == "Quercetin":
                systems=["Neutral"]
            else:
                systems=["Acid", "Neutral"]
            for pH in systems:
                res = open("Rg_results",'w')
                plot = open("plotGnu.plt",'w')

                plot.write("set terminal pngcairo size 2200,600 enhanced font 'Verdana-Bold,10'\n")
                plot.write('\n')
                plot.write('#Setting fonts#\n')
                plot.write('\n')
                plot.write("titleFont=\"'Verdana-Bold,20'\"\n")
                plot.write("labelFont=\"'Verdana-Bold,18'\"\n")
                plot.write("ticsFont=\"'Verdana-Bold,16'\"\n")
                plot.write("keyFont=\"'Verdana-Bold,14'\"\n")
                plot.write('\n')
                plot.write('#Setting the pallete#\n')
                plot.write('\n')
                plot.write('set style line 8 lc rgb "black"         lt 1 dt 1 pt 1 lw 3 ps 1    # error bar / Liu2009 / Wu2010\n')
                plot.write('set style line 1 lc rgb "red"           lt 1 dt 1 pt 1 lw 3 ps 1   # 2016H66\n')
                plot.write('set style line 2 lc rgb "blue"          lt 1 dt 1 pt 1 lw 3 ps 1    # Maingi2012\n')
                plot.write('set style line 3 lc rgb "green"         lt 1 dt 1 pt 1 lw 3 ps 1   # Maiti2004 / Jain2013\n')
                plot.write('set style line 4 lc rgb "dark-green"    lt 1 dt 1 pt 1 lw 3 ps 1   # Maiti2005 \n')
                plot.write('set style line 5 lc rgb "dark-cyan"     lt 1 dt 1 pt 1 lw 3 ps 1    # Opitz2006 \n')
                plot.write('set style line 6 lc rgb "#9a0099"       lt 1 dt 1 pt 1 lw 3 ps 1    # Lee2002   (dark-purple)\n')
                plot.write('set style line 7 lc rgb "#f000f2"       lt 1 dt 1 pt 1 lw 3 ps 1    # Caballero2013 (pink)\n')
                plot.write('set style line 9 lc rgb "orange"        lt 1 dt 1 pt 1 lw 3 ps 1    # Barraza2018\n')
                plot.write('set style line 10 lc rgb "#006b00"      lt 1 dt 1 pt 6 lw 3 ps 2    # Gromos-Kanchi2018\n')
                plot.write('set style line 11 lc rgb "red"          lt 1 dt 1 pt 6 lw 3 ps 2    # CHARMM-Kanchi2018\n')
                plot.write('set style line 12 lc rgb "dark-green"   lt 1 dt 2 pt 8 lw 3 ps 2   # Porcar2008 / Topp1998\n')
                plot.write('set style line 13 lc rgb "blue"         lt 1 dt 2 pt 8 lw 3 ps 2   # Prosa1997\n')
                plot.write('set style line 14 lc rgb "dark-orange"  lt 1 dt 2 pt 8 lw 3 ps 2   # Rathgeber2002 / Scherrenberg1998\n')
                plot.write('set style line 15 lc rgb "#646464"      lt 1 dt 1 pt 1 lw 3 ps 1    # Freire2016 (dark-gray)\n')
                plot.write('set style line 16 lc rgb "dark-red"     lt 1 dt 1 pt 1 lw 3 ps 1    # Tanis2009\n')
                plot.write('\n')
                plot.write('#####################\n')
                plot.write('\n')
                plot.write('set key left font @keyFont \n')
                plot.write('\n')
                plot.write('set border front lc 0 lt 1 lw 4\n')
                plot.write('set grid back linecolor 0 lt 0 lw 2\n')
                plot.write('\n')
                plot.write('#set multiplot layout 1,3\n')
                plot.write('\n')
                plot.write('#-------------------------------//------------------------------#\n')
                plot.write('\n')
                plot.write('set output "RDF.png"\n')
                plot.write('set xlabel "distance(nm)" font @labelFont\n')
                plot.write('set xtics font @ticsFont\n')
                plot.write('#set xrange [0:5]\n')
                plot.write('\n')
                plot.write('set ylabel "g(r)" font @labelFont\n')
                plot.write('set ytics font @ticsFont\n')
                plot.write('#set yrange [0:4]\n')
                plot.write('set title "RDF" font @titleFont\n')
                plot.write('\n')
                plot.write('set output "dist.png"\n')
                plot.write('set xlabel "distance(nm)" font @labelFont\n')
                plot.write('set xtics font @ticsFont\n')
                plot.write('#set xrange [0:5]\n')
                plot.write('\n')
                plot.write('set ylabel "g(r)" font @labelFont\n')
                plot.write('set ytics font @ticsFont\n')
                plot.write('#set yrange [0:4]\n')
                plot.write('set title "# of ligands" font @titleFont\n')
                plot.write('set title "number of ligands within the dendrimer" font @titleFont\n')

                current_case=root+"/{0}/G{1}/{2}/tmp".format(case, G, pH)
                os.chdir(current_case)
                os.system("pwd")

                if pH[0:4] == "Acid":
                    pH="Acid"
                elif pH[0:7] == "Neutral":
                    pH="Neutral"
                else:
                    print("ALGUMA COISA ERRADA\n")
                    exit()

                os.system("mkdir -p {}/proc".format(current_case))
                index(ligand[case])

                Rgm, Rgxm, Rgym, Rgzm = calculate_gyrate(current_case+"/gyrate_{0}G{1}_{2}.xvg".format(case, G, pH), ff, lf)
                res.write("\n\n<------------>{0}G{1}_{2}<------------>\n".format(case, G, pH))
                res.write("\nRg: {0} +/- {1}\n".format(Rgm[0], Rgm[1]))

                # calculate_delta()
                # res.write("\ndelta: {0} +/- {1}\n".format(delta[0], delta[1]))
                
                rdf(current_case+"/rdf_{0}G{1}_{2}.xvg".format(case, G, pH), ff, lf)
                plot.write('"{0}/proc/rdf_{1}G{2}_{3}.xvg" using 1:2 title "{1}G{2}-{3}" with lines ls 0, \\'.format(current_case, case, G, pH))
                
                time, ligands_n = distance(current_case+"/distances_{0}G{1}_{2}.xvg".format(case, G, pH), ff, lf, Rgm[0]+2*Rgm[1])
                file = open("ligands_{0}G{1}_{2}.xvg".format(case, G, pH),'w')
                for t in range(len(time)):
                    file.write("{0:10.3f} {1:10d} \n".format(time[t], ligands_n[t]))
                file.close()
                plot.write('"{0}/proc/ligands_{1}G{2}_{3}.xvg" using 1:2 title "{1}G{2}-{3}" with lines ls 0, \\'.format(current_case, case, G, pH))

                os.system("mv gyrate_{0}G{1}_{2}.xvg rdf_{0}G{1}_{2}.xvg distances_{0}G{1}_{2}.xvg ligands_{0}G{1}_{2}.xvg proc/.".format(case, G, pH))
                clean(["\#*"])

                os.chdir(root)
                
                # plt.plot(time, ligands_n)
                # plt.show()
    res.close()
    plot.close()

if __name__ == "__main__":
    main()
