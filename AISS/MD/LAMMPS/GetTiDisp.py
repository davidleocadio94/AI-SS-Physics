import numpy as np
from scipy import fftpack
from scipy import signal
import matplotlib.pyplot as plt
import glob

def read_dump(FileName,NumAtoms,Nfield, SnapShots):
    """read lammps dumpfile with header (not saved)"""

    File = open(FileName,'r')
    data = np.ndarray((NumAtoms,Nfield,SnapShots),dtype=float)
    
    t = 0
    while (t < SnapShots):
        #read header
        h1 = File.readline()
        time = File.readline()
        h2 = File.readline()
        numatoms = File.readline()
        h3 = File.readline()
        xlen = File.readline()
        ylen = File.readline()
        zlen = File.readline()
        h4 = File.readline()

        for a in range(NumAtoms):
            #Read string -> strip '\n' char -> split into new list
            line = File.readline().strip('\n').split()
            data[a,:,t] = line
            
        t += 1

    File.close()
    return data


def CalculateAngle(v1,v2):
    angle=np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))
    angle=np.clip(angle,-1,1)
    angle=np.arccos(angle)
    return angle

def CalculateDisp(v1,v2):
    delta=np.abs(np.linalg.norm(v1-v2))
    return delta


def GetDispList(num_atoms,NSnapShots,Path):
    #Path = "/home/villarreal/Documents/LAMMPS/135Atoms_RDFAndVautocorr/4x4x4/1585.0/dump.velocity"
    #total = 20000
    #num_atoms = 320
    nfield = 8
    #def read_dump(FileName,NumAtoms,Nfield, SnapShots):
    
    data = read_dump(Path,num_atoms,nfield,NSnapShots)
    
    #initial=data[:,2:5,0]
    for i in range(num_atoms):
        if data[i,0,0]==2:
            v_Ti_i=data[i,2:5,0]
    idxs_init=data[:,0,0]
    disp_list=np.array([])
    k=10000
    snapshots=np.random.randint(1,NSnapShots,k)
    data=data[:,:,snapshots]
    #print(np.shape(data))
    for snapshot in range(k):
     for i in range(num_atoms):
        #print(data[i,1,0])
        #if data[i,1,snapshot]==3 and data[i-1,0,snapshot]==2:#idxs_init[i]:
        if data[i,0,snapshot]==2:
          v_Ti_f=data[i,2:5,snapshot]
          #v_Ti_i=initial[i-1,:]; v_Ti_f=data[i-1,2:5,snapshot]
          #v_O_i=initial[i,:]
          #v_O_f=data[i,2:5,snapshot]
          #v1=v_O_i-v_Ti; v2=v_O_f-v_Ti
          #angle=CalculateAngle(v1,v2)
          disp=CalculateDisp(v_Ti_i,v_Ti_f)

         #print("type",data[i,1,snapshot])
         #print("id",data[i,0,snapshot])
         #print("idxs",idxs_init[i])
         #angle=np.dot(data[i,2:5,snapshot],initial[i,:])/(np.linalg.norm(data[i,2:5,snapshot])*\
         #   np.linalg.norm(initial[i,:]))
         #angle = np.clip(angle, -1, 1)
    
         #angle=np.arccos(angle)
          disp_list=np.append(disp_list,disp)
    
    return disp_list




disp_list_highT=GetDispList(320,20000,"1200/dump.velocity")
disp_list_lowT=GetDispList(320,20000,"140/dump.velocity")


plt.figure()
bins = np.linspace(0, .5, 50)
plt.hist(disp_list_highT,bins,alpha=0.5,label="1200K",density=True)
plt.hist(disp_list_lowT,bins,alpha=0.5,label='140K',density=True)
plt.ylabel("Count")
plt.xlabel("Ti Displacement")
# plt.hist(x, bins, alpha=0.5, label='x')
# plt.hist(y, bins, alpha=0.5, label='y')
plt.legend(loc='upper right')
plt.savefig("1200vs140TiDisp.png")
