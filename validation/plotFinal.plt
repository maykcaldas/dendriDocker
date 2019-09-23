set terminal pngcairo size 1000,1000 enhanced font 'Verdana-Bold,10'

#Setting fonts#

titleFont="'Verdana-Bold,20'"
labelFont="'Verdana-Bold,18'"
ticsFont="'Verdana-Bold,16'"
keyFont="'Verdana-Bold,14'"

#Setting the pallete#

set style line 8 lc rgb "black"         lt 1 dt 1 pt 1 lw 3 ps 1    # error bar / Liu2009 / Wu2010
set style line 1 lc rgb "red"           lt 1 dt 1 pt 1 lw 3 ps 1   # 2016H66
set style line 2 lc rgb "blue"          lt 1 dt 1 pt 1 lw 3 ps 1    # Maingi2012
set style line 3 lc rgb "green"         lt 1 dt 1 pt 1 lw 3 ps 1   # Maiti2004 / Jain2013
set style line 4 lc rgb "dark-green"    lt 1 dt 1 pt 1 lw 3 ps 1   # Maiti2005 
set style line 5 lc rgb "dark-cyan"     lt 1 dt 1 pt 1 lw 3 ps 1    # Opitz2006 
set style line 6 lc rgb "#9a0099"       lt 1 dt 1 pt 1 lw 3 ps 1    # Lee2002   (dark-purple)
set style line 7 lc rgb "#f000f2"       lt 1 dt 1 pt 1 lw 3 ps 1    # Caballero2013 (pink)
set style line 9 lc rgb "orange"        lt 1 dt 1 pt 1 lw 3 ps 1    # Barraza2018
set style line 10 lc rgb "#006b00"      lt 1 dt 1 pt 6 lw 3 ps 2    # Gromos-Kanchi2018
set style line 11 lc rgb "red"          lt 1 dt 1 pt 6 lw 3 ps 2    # CHARMM-Kanchi2018
set style line 12 lc rgb "dark-green"   lt 1 dt 2 pt 8 lw 3 ps 2   # Porcar2008 / Topp1998
set style line 13 lc rgb "blue"         lt 1 dt 2 pt 8 lw 3 ps 2   # Prosa1997
set style line 14 lc rgb "dark-orange"  lt 1 dt 2 pt 8 lw 3 ps 2   # Rathgeber2002 / Scherrenberg1998
set style line 15 lc rgb "#646464"      lt 1 dt 1 pt 1 lw 3 ps 1    # Freire2016 (dark-gray)
set style line 16 lc rgb "dark-red"     lt 1 dt 1 pt 1 lw 3 ps 1    # Tanis2009

#####################

set key left font @keyFont 

set border front lc 0 lt 1 lw 4
set grid back linecolor 0 lt 0 lw 2

#set multifile layout 1,3

#-------------------------------//------------------------------#

set output "RDF.png"
set xlabel "distance(nm)" font @labelFont
set xtics font @ticsFont
#set xrange [0:5]

set ylabel "g(r)" font @labelFont
set ytics font @ticsFont
set yrange [0:50]
plot "/home/mayk/Documents/Labmmol/Dendrimer/dendriDocker/validation/RESULTS/proc/rdfLig_SilybinAG2_n20.xvg" using 1:2 title "Lig-SilybinAG2-n20" with lines ls 1, \
"/home/mayk/Documents/Labmmol/Dendrimer/dendriDocker/validation/RESULTS/proc/rdfWat_SilybinAG2_n20.xvg" using 1:2 title "Water-SilybinAG2-n20" with lines ls 3, \
"/home/mayk/Documents/Labmmol/Dendrimer/dendriDocker/validation/RESULTS/proc/rdfTer_SilybinAG2_n20.xvg" using 1:2 title "Ter-SilybinAG2-n20" with lines ls 2, \
"/home/mayk/Documents/Labmmol/Dendrimer/dendriDocker/validation/RESULTS/proc/rdfDend_SilybinAG2_n20.xvg" using 1:2 title "Dend-SilybinAG2-n20" with lines ls 4, \
#"/home/mayk/Documents/Labmmol/Dendrimer/dendriDocker/validation/RESULTS/proc/rdfLig_SilybinAG2_n25.xvg" using 1:2 title "Lig-SilybinAG2-n25" with lines ls 1, \
#"/home/mayk/Documents/Labmmol/Dendrimer/dendriDocker/validation/RESULTS/proc/rdfTer_SilybinAG2_n25.xvg" using 1:2 title "Water-SilybinAG2-n25" with lines ls 2, \
#"/home/mayk/Documents/Labmmol/Dendrimer/dendriDocker/validation/RESULTS/proc/rdfWat_SilybinAG2_n25.xvg" using 1:2 title "Ter-SilybinAG2-n25" with lines ls 3, \
#"/home/mayk/Documents/Labmmol/Dendrimer/dendriDocker/validation/RESULTS/proc/rdfDend_SilybinAG2_n25.xvg" using 1:2 title "Dend-SilybinAG2-n25" with lines ls 4, \

##################################################
set output "dist.png"
set xlabel "time(ns)" font @labelFont
set xtics font @ticsFont
#set xrange [0:5]

set ylabel "number of ligands" font @labelFont
set ytics font @ticsFont
#set yrange [0:20]
plot "/home/mayk/Documents/Labmmol/Dendrimer/dendriDocker/validation/RESULTS/proc/ligands_SilybinAG2_n20.xvg" using ($1/1000):2 title "SilybinAG2-n20" with lines ls 1, \
#"/home/mayk/Documents/Labmmol/Dendrimer/dendriDocker/validation/RESULTS/proc/ligands_SilybinAG2_n25.xvg" using ($1/1000):2 title "SilybinAG2-n25" with lines ls 2, \

##################################################
set output "ligs.png"
set key right font @keyFont 

set xlabel "time(ns)" font @labelFont
set xtics font @ticsFont
#set xrange [0:5]

set ylabel "distance of ligands" font @labelFont
set ytics font @ticsFont
#set yrange [0:30]

#plot for [n=2:*] "/home/mayk/Documents/Labmmol/Dendrimer/dendriDocker/validation/RESULTS/proc/distances_SilybinAG2_n20.xvg" using ($1/1000):n title sprintf("Lig%d", n-1) with lines lw 3, \
#for [n=2:*] "/home/mayk/Documents/Labmmol/Dendrimer/dendriDocker/validation/RESULTS/proc/distances_SilybinAG2_n25.xvg" using ($1/1000):n title sprintf("Lig%d", n-1) with lines lw 3, \
