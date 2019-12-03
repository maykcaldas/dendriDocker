#!/usr/bin/python3
#!-*- coding: utf8 -*-

'''
Software general documentation

gromacsBuilding was made using python 3.5.2

'''

import os

def write_mdp_em(file_name, emtol, nsteps, emstep):
    if os.path.isfile(file_name):
        print('!!!Backing up the existing em file!!!')
        os.rename(file_name, 'bck.'+file_name)
    em=open(file_name,'w')
    
    #
    em.write('; minim.mdp - used as input into grompp to generate em.tpr\n')
    em.write('integrator	    = steep		; Algorithm (steep = steepest descent minimization)\n')
    em.write('emtol		        = {emtol}  	; Stop minimization when the maximum force < 1000.0 kJ/mol/nm\n'.format(emtol=emtol))
    em.write('emstep      	    = {emstep}      ; Energy step size\n'.format(emstep=emstep))
    em.write('nsteps		    = {nsteps}	  	; Maximum number of (minimization) steps to perform\n'.format(nsteps=nsteps))
    em.write('nstcgsteep	    = 10\n')
    em.write('\n')
    em.write('; Parameters describing how to find the neighbors of each atom and how to calculate the interactions\n')
    em.write('nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces\n')
    em.write('cutoff-scheme     = Verlet\n')
    em.write('ns_type		    = grid		; Method to determine neighbor list (simple, grid)\n')
    em.write('coulombtype	    = PME		; Treatment of long range electrostatic interactions\n')
    em.write('rcoulomb	        = 1.2		; Short-range electrostatic cut-off\n')
    em.write('; VDW\n')
    em.write('vdwtype			= PME\n')
    em.write('fourierspacing    = 0.12		; grid spacing for FFT\n')
    em.write('fourier_nx        = 0\n')
    em.write('fourier_ny        = 0\n')
    em.write('fourier_nz        = 0\n')
    em.write('pme_order         = 4			; cubic interpolation\n')
    em.write('ewald_rtol        = 1e-5\n')
    em.write('rvdw		        = 1.2		; Short-range Van der Waals cut-off\n')
    em.write('pbc		        = xyz 		; Periodic Boundary Conditions (yes/no)\n')
    em.write('rlist             = 1.2\n')


def write_mdp_nvt(file_name, nsteps, dt, cont, temp):
    if os.path.isfile(file_name):
        print('!!!Backing up the existing nvt file!!!')
        os.rename(file_name, 'bck.'+file_name)
    nvt=open(file_name,'w')
    
    nvt.write('title		= PAMAM NVT equilibration\n')
    nvt.write('; Run parameters\n')
    nvt.write('integrator	        = md		; leap-frog integrator\n')
    nvt.write('nsteps		        = {nsteps}	; {time} fs * {nsteps} = {run} ns\n'.format(nsteps=nsteps, time=float(dt)*1000, run=float(nsteps)*float(dt)/1000))
    nvt.write('dt		            = {dt}	    ; {time} fs\n'.format(dt=dt, time=float(dt)*1000))
    nvt.write('; Output control\n')
    nvt.write('nstxout	    	    = 500		; save coordinates every 1.0 ps\n')
    nvt.write('nstvout   		    = 500000    ; save velocities every 1000.0 ps\n')
    nvt.write('nstenergy	        = 500000    ; save energies every 1000.0 ps\n')
    nvt.write('nstlog		        = 500		; update log file every 1.0 ps\n')
    nvt.write('nstxout-compressed   = 500       ; save compressed coordinates every 1.0 ps\n')
    nvt.write('; Bond parameters\n')
    nvt.write('continuation	        = {cont}	; first dynamics run\n'.format(cont=cont))
    nvt.write('constraint_algorithm = lincs	    ; holonomic constraints \n')
    nvt.write('constraints	        = all-bonds	; (even heavy atom-H bonds) constrained\n')
    nvt.write('lincs_iter	        = 1		    ; accuracy of LINCS\n')
    nvt.write('lincs_order	        = 4		    ; also related to accuracy\n')
    nvt.write('; Neighborsearching\n')
    nvt.write('cutoff-schnvte       = Verlet\n')
    nvt.write('rlist                = 1.2\n')
    nvt.write('ns_type		        = grid		; search neighboring grid cells\n')
    nvt.write('nstlist		        = 10	    ; 20 fs, largely irrelevant with Verlet schnvte\n')
    nvt.write('rcoulomb	            = 1.2		; short-range electrostatic cutoff (in nm)\n')
    nvt.write('rvdw		            = 1.2		; short-range van der Waals cutoff (in nm)\n')
    nvt.write('; Electrostatics\n')
    nvt.write('coulombtype	        = PME		; Particle Mesh Ewald for long-range electrostatics\n')
    nvt.write('; VDW\n')
    nvt.write('vdwtype			    = cut-off\n')
    nvt.write('fourierspacing       = 0.12		; grid spacing for FFT\n')
    nvt.write('fourier_nx           = 0\n')
    nvt.write('fourier_ny           = 0\n')
    nvt.write('fourier_nz           = 0\n')
    nvt.write('pme_order            = 4		; cubic interpolation\n')
    nvt.write('ewald_rtol           = 1e-5\n')
    nvt.write('; Tnvtperature coupling is on\n')
    nvt.write('tcoupl		        = V-rescale	            ; modified Berendsen thermostat\n')
    nvt.write('tc-grps		        = non-Water Water		; two coupling groups - more accurate\n')
    nvt.write('tau_t		        = 0.1		0.1		    ; time constant, in ps\n'.format(temp=temp))
    nvt.write('ref_t		        = {temp}	{temp} 		; reference temperature, one for each group, in K\n'.format(temp=temp))
    nvt.write('; Pressure coupling is off\n')
    nvt.write('pcoupl		        = no 		; no pressure coupling in NVT\n')
    nvt.write('; Periodic boundary conditions\n')
    nvt.write('pbc		            = xyz		    ; 3-D PBC\n')
    nvt.write('; Dispersion correction\n')
    nvt.write(';DispCorr	        = EnerPres	; account for cut-off vdW scheme\n')
    nvt.write('; Velocity generation\n')
    if (cont == "no"):
        nvt.write('gen_vel	            = yes		; assign velocities from Maxwell distribution\n')
        nvt.write('gen_temp	            = {temp} ; temperature for Maxwell distribution\n'.format(temp=temp))
        nvt.write('gen_seed	            = -1		; generate a random seed\n')
    elif (cont == "yes"):
        nvt.write('gen_vel	            = no		; assign velocities from Maxwell distribution\n')
        nvt.write('gen_seed	            = -1		; generate a random seed\n')


def write_mdp_ion(file_name):
    if os.path.isfile(file_name):
        print('!!!Backing up the existing ion file!!!')
        os.rename(file_name, 'bck.'+file_name)
    ion=open(file_name,'w')

    ion.write('; ions.mdp - used as input into grompp to generate ions.tpr')
    ion.write('; Parameters describing what to do, when to stop and what to save\n')
    ion.write('integrator	= steep		; Algorithm (steep = steepest descent minimization)\n')
    ion.write('emtol		= 1000.0  	; Stop minimization when the maximum force < 1000.0 kJ/mol/nm\n')
    ion.write('emstep      = 0.01      ; Energy step size\n')
    ion.write('nsteps		= 50000	  	; Maximum number of (minimization) steps to perform\n')
    ion.write('\n')
    ion.write('; Parameters describing how to find the neighbors of each atom and how to calculate the interactions\n')
    ion.write('nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces\n')
    ion.write('cutoff-scheme   = Verlet\n')
    ion.write('ns_type		    = grid		; Method to determine neighbor list (simple, grid)\n')
    ion.write('coulombtype	    = PME		; Treatment of long range electrostatic interactions\n')
    ion.write('rcoulomb	    = 1.0		; Short-range electrostatic cut-off\n')
    ion.write('rvdw		    = 1.0		; Short-range Van der Waals cut-off\n')
    ion.write('pbc		        = xyz 		; Periodic Boundary Conditions (yes/no)\n')


def write_mdp_npt(file_name, nsteps, dt, cont, temp, press):
    if os.path.isfile(file_name):
        print('!!!Backing up the existing npt file!!!')
        os.rename(file_name, 'bck.'+file_name)
    npt=open(file_name,'w')

    npt.write('title		            = PAMAM NPT equilibration \n')
    npt.write('; Run parameters\n')
    npt.write('integrator	            = md		; leap-frog integrator\n')
    npt.write('nsteps		            = {nsteps}	; {time} fs * {nsteps} = {run} ns\n'.format(nsteps=nsteps, time=float(dt)*1000, run=float(dt)*float(nsteps)/1000))
    npt.write('dt		                = {dt}		; {time} fs\n'.format(dt=dt, time=float(dt)*1000))
    npt.write('; Output control\n')
    npt.write('nstxout		            = 500		; save coordinates every 1.0 ps\n')
    npt.write('nstvout		            = 500000	; save velocities every 1000.0 ps\n')
    npt.write('nstenergy	            = 500000	; save energies every 1000.0 ps\n')
    npt.write('nstlog		            = 500		; update log file every 1.0 ps\n')
    npt.write('nstxout-compressed       = 500      ; save compressed coordinates every 1.0 ps\n')
    npt.write('; Bond parameters\n')
    npt.write('continuation	            = {cont}		; Restarting after NVT '.format(cont=cont))
    npt.write('constraint_algorithm     = lincs	    ; holonomic constraints \n')
    npt.write('constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained\n')
    npt.write('lincs_iter	            = 1		    ; accuracy of LINCS\n')
    npt.write('lincs_order	            = 4		    ; also related to accuracy\n')
    npt.write('; Neighborsearching\n')
    npt.write('cutoff-scheme            = Verlet\n')
    npt.write('rlist                    = 1.2\n')
    npt.write('ns_type		            = grid		; search neighboring grid cells\n')
    npt.write('nstlist		            = 10	    ; 20 fs, largely irrelevant with Verlet scheme\n')
    npt.write('rcoulomb	                = 1.2		; short-range electrostatic cutoff (in nm)\n')
    npt.write('rvdw		                = 1.2		; short-range van der Waals cutoff (in nm)\n')
    npt.write('; Electrostatics\n')
    npt.write('coulombtype	            = PME		; Particle Mesh Ewald for long-range electrostatics\n')
    npt.write('; VDW\n')
    npt.write('vdwtype		            = cut-off\n')
    npt.write('fourierspacing           = 0.12		; grid spacing for FFT\n')
    npt.write('fourier_nx               = 0\n')
    npt.write('fourier_ny               = 0\n')
    npt.write('fourier_nz               = 0\n')
    npt.write('pme_order                = 4			; cubic interpolation\n')
    npt.write('ewald_rtol               = 1e-5\n')
    npt.write('; Temperature coupling is on\n')
    npt.write('tcoupl		            = V-rescale	            ; modified Berendsen thermostat\n')
    npt.write('tc-grps		            = non-Water Water		; two coupling groups - more accurate\n')
    npt.write('tau_t		            = 0.1	  	0.1		    ; time constant, in ps\n')
    npt.write('ref_t		            = {temp}       {temp} 		; reference temperature, one for each group, in K'.format(temp=temp))
    npt.write('; Pressure coupling is on\n')
    npt.write('pcoupl		            = Berendsen	    ; Pressure coupling on in NPT\n')
    npt.write('pcoupltype	            = isotropic	            ; uniform scaling of box vectors\n')
    npt.write('tau_p		            = 2.0		            ; time constant, in ps\n')
    npt.write('ref_p		            = {press}		            ; reference pressure, in bar'.format(press=press))
    npt.write('compressibility          = 4.5e-5	            ; isothermal compressibility of water, bar^-1\n')
    npt.write('refcoord_scaling         = com\n')
    npt.write('; Periodic boundary conditions\n')
    npt.write('pbc		                = xyz		; 3-D PBC\n')
    npt.write('; Dispersion correction\n')
    npt.write(';DispCorr	            = EnerPres	; account for cut-off vdW scheme\n')
    npt.write('; Velocity generation\n')
    if (cont == "no"):
        npt.write('gen_vel	                = yes		; assign velocities from Maxwell distribution\n')
        npt.write('gen_temp	                = {temp} ; temperature for Maxwell distribution\n'.format(temp=temp))
        npt.write('gen_seed	                = -1		; generate a random seed\n')
    elif (cont == "yes"):
        npt.write('gen_vel	                = no		; assign velocities from Maxwell distribution\n')
        npt.write('gen_seed	                = -1		; generate a random seed\n')


def write_mdp_md(file_name, nsteps, dt, temp, press):
    if os.path.isfile(file_name):
        print('!!!Backing up the existing md file!!!')
        os.rename(file_name, 'bck.'+file_name)
    md=open(file_name,'w')

    md.write('title		= PAMAM Dendrimer G0 in Basic environment simulation\n') 
    md.write('; Run parameters\n')
    md.write('integrator	            = md		; leap-frog integrator\n')
    md.write('nsteps		            = {nsteps}	; {time} fs * {nsteps} = {run} ns\n'.format(nsteps=nsteps, time=float(dt)*1000, run=float(dt)*float(nsteps)/1000))
    md.write('dt		                = {dt}		; {time} fs\n'.format(dt=dt, time=float(dt)*1000))
    md.write('; Output control\n')
    md.write('nstxout		            = 5000		; save coordinates every 10.0 ps\n')
    md.write('nstvout		            = 500000	; save velocities every 1000.0 ps\n')
    md.write('nstenergy	                = 500000	; save energies every 1000.0 ps\n')
    md.write('nstlog		            = 500000	; update log file every 1000.0 ps\n')
    md.write('nstxout-compressed        = 5000      ; save compressed coordinates every 10.0 ps\n')
    md.write('                                ; nstxout-compressed replaces nstxtcout\n')
    md.write('compressed-x-grps         = System    ; replaces xtc-grps\n')
    md.write('; Bond parameters\n')
    md.write('continuation	            = yes		; Restarting after NPT \n')
    md.write('constraint_algorithm      = lincs	    ; holonomic constraints \n')
    md.write('constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained\n')
    md.write('lincs_iter	            = 1		    ; accuracy of LINCS\n')
    md.write('lincs_order	            = 2		    ; also related to accuracy\n')
    md.write('; Neighborsearching\n')
    md.write('cutoff-scheme             = Verlet\n')
    md.write('rlist                     = 1.2\n')
    md.write('ns_type		            = grid		; search neighboring grid cells\n')
    md.write('nstlist		            = 10	    ; 20 fs, largely irrelevant with Verlet scheme\n')
    md.write('rcoulomb	                = 1.2		; short-range electrostatic cutoff (in nm)\n')
    md.write('rvdw		                = 1.2		; short-range van der Waals cutoff (in nm)\n')
    md.write('; Electrostatics\n')
    md.write('coulombtype	            = PME		; Particle Mesh Ewald for long-range electrostatics\n')
    md.write('; VDW\n')
    md.write('vdwtype	                = cut-off\n')
    md.write('fourierspacing            = 0.12		; grid spacing for FFT\n')
    md.write('fourier_nx                = 0\n')
    md.write('fourier_ny                = 0\n')
    md.write('fourier_nz                = 0\n')
    md.write('pme_order                 = 4			; cubic interpolation\n')
    md.write('ewald_rtol                = 1e-5\n')
    md.write('; Temperature coupling is on\n')
    md.write('tcoupl	                = V-rescale	            ; modified Berendsen thermostat\n')
    md.write('tc-grps	                = non-Water Water		; two coupling groups - more accurate\n')
    md.write('tau_t		                = 0.1		0.1	        ; time constant, in ps\n')
    md.write('ref_t		                = {temp}	{temp}      ; reference temperature, one for each group, in K\n'.format(temp=temp))
    md.write('; Pressure coupling is on\n')
    md.write('pcoupl		            = Parrinello-Rahman	    ; Pressure coupling on in NPT\n')
    md.write('pcoupltype	            = isotropic	            ; uniform scaling of box vectors\n')
    md.write('tau_p		                = 2.0		            ; time constant, in ps\n')
    md.write('ref_p		                = {press}		            ; reference pressure, in bar\n'.format(press=press))
    md.write('compressibility           = 4.5e-5	            ; isothermal compressibility of water, bar^-1\n')
    md.write('; Periodic boundary conditions\n')
    md.write('pbc		                = xyz		; 3-D PBC\n')
    md.write('; Dispersion correction\n')
    md.write(';DispCorr	                = EnerPres	; account for cut-off vdW scheme\n')
    md.write('; Velocity generation\n')
    md.write('gen_vel		            = no		; Velocity generation is off \n')
    md.write('gen_seed	                = -1		; generate a random seed\n')


def write_mdp(workflow):
    for step in workflow:
        if step[0] == "BOX":
            pass
        elif step[0] == "SOLV":
            pass
        elif step[0] == "EM":
            em_file_name=step[1]["file_name"]
            em_emtol=step[1]["emtol"]
            em_nsteps=step[1]["nsteps"]
            em_emstep=step[1]["emstep"]
            
            write_mdp_em(em_file_name, em_emtol, em_nsteps, em_emstep)
        elif step[0] == "ION":
            ion_file_name=step[1]["file_name"]

            write_mdp_ion(ion_file_name)
        elif step[0] == "INSERT":
            pass
        elif step[0] == "SD":
            write_mdp_sd()
        elif step[0] == "NVT":
            nvt_file_name=step[1]["file_name"]
            nvt_dt=step[1]["dt"]
            nvt_nsteps=step[1]["nsteps"]
            nvt_cont=step[1]["cont"]
            nvt_temp=step[1]["temp"]
            
            write_mdp_nvt(nvt_file_name, nvt_nsteps, nvt_dt, nvt_cont, nvt_temp)
        elif step[0] == "NPT":
            npt_file_name=step[1]["file_name"]
            npt_dt=step[1]["dt"]
            npt_nsteps=step[1]["nsteps"]
            npt_cont=step[1]["cont"]
            npt_temp=step[1]["temp"]
            npt_press=step[1]["press"]
            
            write_mdp_npt(npt_file_name, npt_nsteps, npt_dt, npt_cont, npt_temp, npt_press)
        elif step[0] == "MD":
            md_file_name=step[1]["file_name"]
            md_dt=step[1]["dt"]
            md_nsteps=step[1]["nsteps"]
            md_cont=step[1]["cont"]
            md_temp=step[1]["temp"]
            md_press=step[1]["press"]
            
            write_mdp_md(md_file_name, md_nsteps, md_dt, md_temp, md_press)
        else:
            error("There is not such mdp process implemented. Please, check your input or contact the developers.")

def error(message):
    print("An error occurred.")
    print(message)
    print("See './__main__.py -h' for more instructions.")
    exit()