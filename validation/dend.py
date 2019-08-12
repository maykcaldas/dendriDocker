import math
import sys
import os
import numpy

f=str(sys.argv[1])
rmax=0.0 #nm
nmin=float(sys.argv[2]) #40000.0 #sys.argv[2]
nmax=float(sys.argv[3]) #50000.0 #sys.argv[3]

def readFile(f):
    try:
        file=open(f,"r")
    except IOError:
        print ("This shit \'" + str(f) +"\' isn't opening")
        exit()
        
    gro=file.readlines()
    file.close()

    #for i in gro: print(i[:-1])

    return(gro)


def readFrame(gro, frame, natoms):

    atoms=[]
    
    firstLine=(int(frame)*(int(natoms)+3))
    lastLine=(int(frame)*(int(natoms)+3)+(int(natoms)+2))

    thisFrame=gro[firstLine:lastLine+1]

    title=(gro[firstLine])
    read=gro[(firstLine+2):(lastLine)]
    box=gro[lastLine]

    for i in range(len(read)):
        atoms.append(read[i].split())
    return(atoms, title)
    

def selDend(atoms):
    groups=['CORE','INTR','TER']
    dend=[]
    for atom in atoms: 
        #print (atom[0].find("CORE"))
        #print ([g in atom[0] for g in groups])
        if (any([g in atom[0] for g in groups]) != False):
            dend.append(atom)
            
    return(dend)

def selMass(atom):
    CMass=14
    NMass=14
    OMass=16
    HMass=1
    
    if atom[1][0] == 'C':
        m = CMass
    elif atom[1][0] == 'N':
        m = NMass
    elif atom[1][0] == 'O':
        m = OMass
    elif atom[1][0] == 'H':
        m = HMass
    else:
        print (atom)
        print ("This is not a valid atom (sorry, we're lazy. Our FF isn't that big')")
        exit()
        
    return (m)


def calcCM(group):
    m=[]
    x=[]
    y=[]
    z=[]
    CMx=0
    CMy=0
    CMz=0
    TMass=0
    for atom in group: # I should have done the sum on this first loop
        m.append(selMass(atom))
        x.append(float(atom[3]))
        y.append(float(atom[4]))
        z.append(float(atom[5]))
        
    for k in range((len(group))):
        TMass+=m[k]
        CMx+=m[k]*x[k]
        CMy+=m[k]*y[k]
        CMz+=m[k]*z[k]
    
    #CM=math.sqrt(CMx*CMx+CMy*CMy+CMz*CMz)
        
    return([CMx/TMass, CMy/TMass, CMz/TMass])


def calcRg(dend, time):
    Rgx2=0
    Rgy2=0
    Rgz2=0
    TMass=0
    
    CM=calcCM(dend)
    for atom in dend:
        mass = selMass(atom)
        Rgx2+=(mass*(float(atom[3])-CM[0])*(float(atom[3])-CM[0]))
        Rgy2+=(mass*(float(atom[4])-CM[1])*(float(atom[4])-CM[1]))
        Rgz2+=(mass*(float(atom[5])-CM[2])*(float(atom[5])-CM[2]))
        TMass+=mass
        
        
    Rg=math.sqrt(Rgx2/TMass+Rgy2/TMass+Rgz2/TMass)
    Rgx=math.sqrt(Rgx2/TMass)
    Rgy=math.sqrt(Rgy2/TMass)
    Rgz=math.sqrt(Rgz2/TMass)
    
    return ([time, Rg, Rgx, Rgy, Rgz])

def writeRg(Rg):
    if os.path.isfile('Rg.txt'):
         os.rename('Rg.txt', 'Rg_bkp.txt')
         fRg=open('Rg.txt','a')
    else:
        fRg=open('Rg.txt','a')
    
    for k in range(9):
        fRg.write("#\n")
    fRg.write("It'll be better when finished. I promisse.\n")
    fRg.write("#\n")
    fRg.write('{0:8} {1:8} {2:8} {3:8} {4:8} \n'.format('time/ps','Rg/nm','Rgx/nm','Rgy/nm','Rgz/nm'))
    fRg.write('@    title "Moments of inertia (total and around axes)"\n')
    fRg.write('@    xaxis  label "Time (ps)"\n')
    fRg.write('@    yaxis  label "I (a.m.u. nm\S2\N)"\n')
    fRg.write('@TYPE xy\n')
    fRg.write('@ view 0.15, 0.15, 0.75, 0.85\n')
    fRg.write('@ legend on\n')
    fRg.write('@ legend box on\n')
    fRg.write('@ legend loctype view\n')
    fRg.write('@ legend 0.78, 0.8\n')
    fRg.write('@ legend length 2\n')
    fRg.write('@ s0 legend "Itot"\n')
    fRg.write('@ s1 legend "I1"\n')
    fRg.write('@ s2 legend "I2"\n')
    fRg.write('@ s3 legend "I3"\n')

    for k in Rg:
        fRg.write('{0:6.2f} {1:8.6f} {2:8.6f} {3:8.6f} {4:8.6f} \n'.format(k[0],k[1],k[2],k[3],k[4]))
        
    fRg.close()


def selSolv(atoms):
    groups=['SOL', 'CL']
    solv=[]
    for atom in atoms: 
        #print (atom[0].find("SOL"))
        #print ([g in atom[0] for g in groups])
        if (any([g in atom[0] for g in groups]) != False):
            solv.append(atom)
            
    return(solv)


def selectMoi(dend, solv, rmax):
    select=[]
    for atom in dend:
        select.append(atom[2])
        for sol in solv:
            dx=float(atom[-3])-float(sol[-3])
            dy=float(atom[-2])-float(sol[-2])
            dz=float(atom[-1])-float(sol[-1])
            r=math.sqrt( (dx)**2 + (dy)**2 + (dz)**2 )
            
            if (sol[2] not in select) and (r<=rmax):
                select.append(sol[2]) #select.append(sol)
                
    return(select)


def calcMoi(select, time):
    G=[[0,0,0],\
       [0,0,0],\
       [0,0,0]]
    TMass=0

    CM=calcCM(select)
    for atom in select:
        mass = (selMass(atom))
        G[0][0]+=(mass*(float(atom[3])-CM[0])*(float(atom[3])-CM[0]))
        G[0][1]+=(mass*(float(atom[3])-CM[0])*(float(atom[4])-CM[1]))
        G[0][2]+=(mass*(float(atom[3])-CM[0])*(float(atom[5])-CM[2]))
        G[1][0]+=(mass*(float(atom[4])-CM[1])*(float(atom[3])-CM[0]))
        G[1][1]+=(mass*(float(atom[4])-CM[1])*(float(atom[4])-CM[1]))
        G[1][2]+=(mass*(float(atom[4])-CM[1])*(float(atom[5])-CM[2]))
        G[2][0]+=(mass*(float(atom[5])-CM[2])*(float(atom[3])-CM[0]))
        G[2][1]+=(mass*(float(atom[5])-CM[2])*(float(atom[4])-CM[1]))
        G[2][2]+=(mass*(float(atom[5])-CM[2])*(float(atom[5])-CM[2]))
        TMass+=mass
    
    for i in G:
        for j in i:
            j=j/TMass
            
    try:
        eigen=numpy.linalg.eig(G)[0]
    except LinAlgError:
        print("The shape tensor eigenvalues calculation didn't converge")
        exit()
        
    eigen.sort()
    
    I=math.sqrt(eigen[0]*eigen[0]+eigen[1]*eigen[1]+eigen[2]*eigen[2])
    
    return ([time, I, eigen[0], eigen[1], eigen[2]])

def writeMoi(moi):
    if os.path.isfile('Moi.txt'):
         os.rename('Moi.txt', 'Moi_bkp.txt')
         fMoi=open('Moi.txt','a')
    else:
        fMoi=open('Moi.txt','a')
        
    for k in range(9):
        fMoi.write("#\n")
    fMoi.write("It'll be better when finished. I promisse.\n")
    fMoi.write("#\n")
    fMoi.write('{0:8} {1:8} {2:8} {3:8} \n'.format('time/ps','Ix/nm2','Iy/nm2','Iz/nm2'))
    fMoi.write('@    title "Moments of inertia (total and around axes)"\n')
    fMoi.write('@    xaxis  label "Time (ps)"\n')
    fMoi.write('@    yaxis  label "I (a.m.u. nm\S2\N)"\n')
    fMoi.write('@TYPE xy\n')
    fMoi.write('@ view 0.15, 0.15, 0.75, 0.85\n')
    fMoi.write('@ legend on\n')
    fMoi.write('@ legend box on\n')
    fMoi.write('@ legend loctype view\n')
    fMoi.write('@ legend 0.78, 0.8\n')
    fMoi.write('@ legend length 2\n')
    fMoi.write('@ s0 legend "Itot"\n')
    fMoi.write('@ s1 legend "I1"\n')
    fMoi.write('@ s2 legend "I2"\n')
    fMoi.write('@ s3 legend "I3"\n')

    for k in moi:
        fMoi.write('{0:6.2f} {1:8.4f} {2:8.4f} {3:8.4f} {4:8.4f}\n'.format(k[0],k[1],k[2],k[3],k[4]))
        
    fMoi.close()


#def calcMean():
#    continue

#def fprintf():
#    continue

def main():
    gro=readFile(f)
    #for i in gro: print (i[0:15])
    natoms=gro[1][:-1]
    #print("natoms:", natoms)
    
    totalFrames = len(gro)/(int(gro[1])+3)
    usedFrames=0
#    print ("Frames Totais:", totalFrames)


    Rg=[]    
    moi=[]
    moix=0
    moiy=0
    moiz=0
    Rgm=0
    for frame in range(totalFrames):
#        print ("Frame:", frame  )
        atoms, title=readFrame(gro, frame, natoms)
        
        time=float(title.split()[-1])

        if (time < nmin):
#            print (str(time) + " Frame was not used")
            continue
        elif (time > nmax):
#            print (str(time) + " Last frame reached")
            break
        else:
            usedFrames+=1
#            print (str(time) + " Computing...")
        
        dend=selDend(atoms)
#        for d in dend: print (d)
        
#        solv=selSolv(atoms)
#        for s in solv: print(s)
        
#        select=selectMoi(dend, solv, rmax)
#        for sel in select: print(sel) #select is a list of indexes

        CM=calcCM(dend)
#        print(CM)
        
        
        asd=calcRg(dend,time)
        Rg.append(calcRg(dend,time))
        Rgm+=asd[1]
#        print(Rg)
        
        moment=calcMoi(dend,time)
        moi.append(moment)
        moix+=moment[-3]
        moiy+=moment[-2]
        moiz+=moment[-1]
        
    writeRg(Rg)
    writeMoi(moi)
    
    Ix=moix/usedFrames
    Iy=moiy/usedFrames
    Iz=moiz/usedFrames
    Rgm=Rgm/usedFrames

    print("Average values:")
    print(Ix, Iy, Iz)
   
    I1=Ix+Iy+Iz
    I2=Ix*Iy+Ix*Iz+Iy*Iz
    
    delta=1.0-(3.0*(I2/(I1*I1)))

    print("Rg: ", math.sqrt(Ix+Iy+Iz), Rgm)
    print("Asphericity: ", delta)
    print("\n")

main()
