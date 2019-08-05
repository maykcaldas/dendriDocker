#!/bin/bash -l 

gmx514=/home/mayk/programs/gromacs-5.1.4/install/bin/gmx
gmx407=/home/mayk/programs/gromacs-4.0.7/install/bin/g_rdf

# 0: Não executa o processamento
# 1: Pós-processa
gyrate=0
moi=0
sasa=0
rdfc=0
rdfd=0
comd=0
dihedral=0
rmsd=0
cutrmsd=0

gro=1

#frames to analysis
ff=10000
lf=100000

#Diffusivity coefficient fit
bgfit=3000
edfit=5000

dend='PAMAM'
setup='CUT'

if [ $cutrmsd == 1 ]; then

    rm diffusionG* final.txt analyze* MSD_*

fi



for i in 0 1 2 3 4 5
do

cd G${i}

    for j in Acid Neutral Basic
    do

        cd ${j}/md

        if [ "${j}" != "Basic" ]; then
            if [ ${i} == 0 ]; then
                $gmx514 make_ndx -f md.gro -o index.ndx << !
                    keep 0
                    r CORE
                    name 1 CORE
                    r TER
                    name 2 TER
                    r SOL
                    name 3 Water
                    r CL-
                    name 4 CL
                    1|2
                    name 5 Dendrimer
                    
                    q
!
            elif [ ${i} != 0 ]; then
                $gmx514 make_ndx -f md.gro -o index.ndx << !
                    keep 0
                    r CORE
                    name 1 CORE
                    r INTR
                    name 2 INTR
                    r TER
                    name 3 TER
                    r SOL
                    name 4 Water
                    r CL-
                    name 5 CL
                    1|2|3
                    name 6 Dendrimer
                    
                    q
!
			fi
        elif [ ${j} == "Basic" ]; then
            if [ ${i} == 0 ]; then
                $gmx514 make_ndx -f md.gro -o index.ndx << !
                    keep 0
                    r CORE
                    name 1 CORE
                    r TER
                    name 2 TER
                    r SOL
                    name 3 Water
                    1|2
                    name 4 Dendrimer
                    
                    q
!
            elif [ ${i} != 0 ]; then
                $gmx514 make_ndx -f md.gro -o index.ndx << !
                    keep 0
                    r CORE
                    name 1 CORE
                    r INTR
                    name 2 INTR
                    r TER
                    name 3 TER
                    r SOL
                    name 4 Water
                    1|2|3
                    name 5 Dendrimer
                    
                    q
!
            fi
        fi
                
		if [ $gyrate == 1 ]; then
			echo Dendrimer | $gmx514 gyrate -s md.tpr \
											-f md.trr \
											-n index.ndx \
											-b ${ff} \
											-e ${lf} \
											-o gyrate_G${i}_${j}.xvg
			mv gyrate_G${i}_${j}.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
		fi
		
		if [ $moi == 1 ]; then
			echo Dendrimer | $gmx514 gyrate -moi \
											-s md.tpr \
											-f md.trr \
											-n index.ndx \
											-b ${ff} \
											-e ${lf} \
											-o moi_G${i}_${j}.xvg
			mv moi_G${i}_${j}.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
		fi
		
		if [ $sasa == 1 ]; then
			echo Dendrimer | $gmx514 sasa  	-s md.tpr \
											-f md.trr \
											-n index.ndx \
											-b ${ff} \
											-e ${lf} \
											-probe 0.14 \
											-o sasa_G${i}_${j}.xvg
			mv sasa_G${i}_${j}.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
		fi
		
		if [ $gro == 1 ]; then
		    echo Dendrimer | $gmx514 trjconv -f md.xtxc -s md.tpr -n index.ndx -o ${dend}_G${i}_${j}.gro -pbc mol 
		fi
		
		mv ${dend}_G${i}_${j}.gro /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.


		if [ $rdfc == 1 ]; then
			echo Dendrimer TER | $gmx407 -f ${dend}_G${i}_${j}.gro \
		 							     -s ${dend}_G${i}_${j}.gro \
		 						    	 -o ter_G${i}_${j}.xvg \
		 					    		 -com \
		 				    			 -n index.ndx \
		 			    				 -b ${ff} \
    	 								 -e ${lf} \
	     								 -bin 0.04
			mv ter_G${i}_${j}.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
			
			if [${j} != "Basic" ]; then
				echo Dendrimer CL | $gmx407 	-f ${dend}_G${i}_${j}.gro \
										-s ${dend}_G${i}_${j}.gro \
										-o cl_G${i}_${j}.xvg \
										-com \
										-n index.ndx \
										-b ${ff} \
										-e ${lf} \
										-bin 0.04
				mv cl_G${i}_${j}.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
			fi

			echo Dendrimer Water | $gmx407 	-f ${dend}_G${i}_${j}.gro \
										-s ${dend}_G${i}_${j}.gro \
										-o water_G${i}_${j}.xvg \
										-com \
										-n index.ndx \
										-b ${ff} \
										-e ${lf} \
										-bin 0.04
			mv water_G${i}_${j}.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
			
			echo Dendrimer Dendrimer | $gmx407 	-f ${dend}_G${i}_${j}.gro \
											-s ${dend}_G${i}_${j}.gro \
											-o rdf_G${i}_${j}.xvg \
											-com \
											-n index.ndx \
											-b ${ff} \
											-e ${lf} \
											-bin 0.04
			mv rdf_G${i}_${j}.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
		fi
		
		if [ $rdfd == 1 ]; then
			echo Dendrimer TER | $gmx407 	-f ${dend}_G${i}_${j}.gro \
											-s ${dend}_G${i}_${j}.gro \
											-o ter_G${i}_${j}_nonorm.xvg \
											-nonorm \
											-com \
											-n index.ndx \
											-b ${ff} \
											-e ${lf} \
											-bin 0.04
			mv ter_G${i}_${j}_nonorm.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
			
			if [ ${j} != "Basic" ]; then
				echo Dendrimer CL | $gmx407 -f ${dend}_G${i}_${j}.gro \
											-s ${dend}_G${i}_${j}.gro \
											-o cl_G${i}_${j}_nonorm.xvg \
											-nonorm \
											-com \
											-n index.ndx \
											-b ${ff} \
											-e ${lf} \
											-bin 0.04
				mv cl_G${i}_${j}_nonorm.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
			fi
		
			echo Dendrimer Water | $gmx407 	-f ${dend}_G${i}_${j}.gro \
											-s ${dend}_G${i}_${j}.gro \
											-o water_G${i}_${j}_nonorm.xvg \
											-nonorm \
											-com \
											-n index.ndx \
											-b ${ff} \
											-e ${lf} \
											-bin 0.04
			mv water_G${i}_${j}_nonorm.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
			
			echo Dendrimer Dendrimer | $gmx407 	-f ${dend}_G${i}_${j}.gro \
												-s ${dend}_G${i}_${j}.gro \
												-o rdf_G${i}_${j}_nonorm.xvg \
												-nonorm \
												-com \
												-n index.ndx \
												-b ${ff} \
												-e ${lf} \
												-bin 0.04
			mv rdf_G${i}_${j}_nonorm.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
		fi
		
		if [ $comd == 1 ]; then
			$gmx514 distance	-select 'com of group CORE plus com of group Dendrimer' \
								-s md.tpr \
								-f md.trr \
								-n index.ndx \
								-oav comdav_G${i}_${j}.xvg \
								-oh comdh_G${i}_${j}.xvg
			mv comdav_G${i}_${j}.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
			mv comdh_G${i}_${j}.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
		fi


		if [ $dihedral == 1 ]; then
    		ltot=$(grep "gd_14" ../${dend}_G${i}_${j}.itp | wc -l) #Counting the number of terminal and internal amides
    		ext=$((2**(${i}+2)))
    		int=$(($ltot-$ext))

			$gmx514 mk_angndx -s md.tpr -n angle.ndx
			echo [ dihedralint ] >> angle.ndx
			grep "gd_14" ../${dend}_G${i}_${j}.itp | tail -$int | awk '{print $1, $2, $3, $4}' >> angle.ndx
			echo "" >> angle.ndx

			echo [ dihedralext ] >> angle.ndx
			grep "gd_14" ../${dend}_G${i}_${j}.itp | tail -$ext | awk '{print $1, $2, $3, $4}' >> angle.ndx
			echo "" >> angle.ndx
			
			#$gmx514 make_ndx -f md.gro -o angle.ndx
			$gmx514 gangle	-g1 dihedral \
							-group1 'group dihedralint' \
							-f ${dend}_G${i}_${j}.gro \
							-s ${dend}_G${i}_${j}.gro \
							-n angle.ndx \
							-binw 2 \
							-oav gangleavint_G${i}_${j}.xvg \
							-oh ganglehint_G${i}_${j}.xvg 
			mv gangleavint_G${i}_${j}.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
			mv ganglehint_G${i}_${j}.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
			
			$gmx514 gangle	-g1 dihedral \
							-group1 'group dihedralext' \
							-f ${dend}_G${i}_${j}.gro \
							-s ${dend}_G${i}_${j}.gro \
							-n angle.ndx \
							-binw 2 \
							-oav gangleavext_G${i}_${j}.xvg \
							-oh ganglehext_G${i}_${j}.xvg 
			mv gangleavext_G${i}_${j}.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
			mv ganglehext_G${i}_${j}.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
		fi

		if [ $rmsd == 1 ]; then
        
        echo Calculating D of the PAMAM_G${i}_${j}

        #echo Dendrimer | $gmx514 msd -f md.xtc\
        #                             -s md.tpr\
        #                             -n index.ndx\
        #                             -b 10000\
        #                             -e 20000\
        #                             -beginfit $bgfit\
        #                             -endfit $edfit\
        #                             -o msd_G${i}_${j}.1.xvg

        #echo Dendrimer | $gmx514 msd -f md.xtc\
        #                             -s md.tpr\
        #                             -n index.ndx\
        #                             -b 20000\
        #                             -e 30000\
        #                             -beginfit $bgfit\
        #                             -endfit $edfit\
        #                             -o msd_G${i}_${j}.2.xvg
        
        #echo Dendrimer | $gmx514 msd -f md.xtc\
        #                             -s md.tpr\
        #                             -n index.ndx\
        #                             -b 30000\
        #                             -e 40000\
        #                             -beginfit $bgfit\
        #                             -endfit $edfit\
        #                             -o msd_G${i}_${j}.3.xvg
        
        #echo Dendrimer | $gmx514 msd -f md.xtc\
        #                             -s md.tpr\
        #                             -n index.ndx\
        #                             -b 40000\
        #                             -e 50000\
        #                             -beginfit $bgfit\
        #                             -endfit $edfit\
        #                             -o msd_G${i}_${j}.4.xvg
        
        echo Dendrimer | $gmx514 msd -f md.xtc\
                                     -s md.tpr\
                                     -n index.ndx\
                                     -b ${ff}\
                                     -e ${lf}\
                                     -beginfit $bgfit\
                                     -endfit $edfit\
                                     -o msd_G${i}_${j}.xvg

        #mv msd_G${i}_${j}.1.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
        #mv msd_G${i}_${j}.2.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
        #mv msd_G${i}_${j}.3.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
        #mv msd_G${i}_${j}.4.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.
        mv msd_G${i}_${j}.xvg /home/mayk/Documents/Labmmol/Dendrimer/Results${dend}_${setup}/proc/.

        fi

		if [ $cutrmsd == 1 ]; then

            cd ../../..
            
            grep "PAMAM\|cm\^2/s" MSD.log > MSD_G${i}_${j}.log

            #grep -A 4 "PAMAM_G${i}_${j}" MSD2.log > MSD_G${i}_${j}.log

            grep -A 1 "PAMAM_G${i}_${j}" MSD_G${i}_${j}.log | cut -c 15-20 | tail -1 > analyze_G${i}_${j}.xvg

            $gmx514 analyze -f analyze_G${i}_${j}.xvg > analyzed_G${i}_${j}.xvg

            echo PAMAM_G${i}_${j} >> final.txt

            grep "SS1" analyzed_G${i}_${j}.xvg >> final.txt

            echo ${i} `grep -A 1 "PAMAM_G${i}_${j}" final.txt | cut -c 5-35 | tail -1 | awk '{print $1*10, $2*10}'` >> diffusionG${j}2.txt

            cd -
        fi

		cd ../..
    done
cd ..
done

./clean.sh
