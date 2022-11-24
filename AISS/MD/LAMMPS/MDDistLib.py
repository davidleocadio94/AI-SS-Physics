import numpy as np
# from scipy import fftpack
# from scipy import signal
from scipy.stats import kurtosis
from scipy.stats import skew
import matplotlib.pyplot as plt
import glob
#from IPython.core.interactiveshell import InteractiveShell
#InteractiveShell.ast_node_interactivity = "all"

class MDFile():


    def __init__(self,PathMD,PathNotTilted,NumAtoms,TotalSnapShotsForAvg,TotalMDLength=30e3):
        
        
        
        self.PathMD=PathMD; self.PathNotTilted=PathNotTilted
        self.NumAtoms=NumAtoms
        self.TotalMDLength=TotalMDLength; self.TotalSnapShotsForAvg=TotalSnapShotsForAvg










    def compute_angle(self,v1,v2,axis):
        v1=v1/np.linalg.norm(v1); v2=v2/np.linalg.norm(v2)
        ex=np.array([1,0,0])
        ey=np.array([0,1,0])
        ez=np.array([0,0,1])
        if axis=="x":
            v1y=np.dot(v1,ey);v2y=np.dot(v2,ey)
            v1z=np.dot(v1,ez);v2z=np.dot(v2,ez)
            v1=np.array([v1y,v1z]);v2=np.array([0,-1])
            v1=v1/np.linalg.norm(v1); v2=v2/np.linalg.norm(v2)
            angle=np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))
            angle=np.clip(angle,-1,1)
            angle=np.arccos(angle)*180/np.pi
            sign=np.sign(np.cross(v1,v2))
            return angle,sign
        elif axis=="y":
            v1x=np.dot(v1,ex);v2x=np.dot(v2,ex)
            v1z=np.dot(v1,ez);v2z=np.dot(v2,ez)
            v1=np.array([v1x,v1z]);v2=np.array([-1,0])
            v1=v1/np.linalg.norm(v1); v2=v2/np.linalg.norm(v2)
            angle=np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))
            angle=np.clip(angle,-1,1)
            angle=np.arccos(angle)*180/np.pi
            sign=np.sign(np.cross(v1,v2))
            return angle,sign
        elif axis=="z":
            v1x=np.dot(v1,ex);v2x=np.dot(v2,ex)
            v1y=np.dot(v1,ey);v2y=np.dot(v2,ey)
            v1=np.array([v1x,v1y]);v2=np.array([0,-1])
            v1=v1/np.linalg.norm(v1); v2=v2/np.linalg.norm(v2)
            angle=np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))
            angle=np.clip(angle,-1,1)
            angle=np.arccos(angle)*180/np.pi
            sign=np.sign(np.cross(v1,v2))
            return angle,sign






    def read_dump(self,FileName,NumAtoms,Nfield, SnapShots):
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



    def GetAngleDistMD(self):

        dataMD = self.read_dump(self.PathMD,self.NumAtoms,8,self.TotalMDLength)
        #dataNotTilted=self.read_dump(self.PathNotTilted,self.NumAtoms,5,1)
        dataNotTilted=self.read_dump(self.PathMD,self.NumAtoms,8,1)
        #dataNotTilted=dataMD[:,:,0]
        angle_list_x=np.array([]);angle_list_y=np.array([]);angle_list_z=np.array([])
        d_list_x=np.array([]);d_list_y=np.array([]);d_list_z=np.array([])
        #snapshots=np.random.randint(1,self.TotalSnapShotsForAvg,self.TotalMDLength)
        snapshots=np.random.randint(1,self.TotalMDLength,self.TotalSnapShotsForAvg)
        dataMD=dataMD[:,:,snapshots]
        #dataMD=np.sort(dataMD,axis=1)
        # for i in np.arange(np.shape(dataMD)[0]): print(dataMD[i,:,0])
        #dataNotTilted=np.sort(dataNotTilted,axis=0)

        dataNotTilted=dataNotTilted[dataNotTilted[:,0,0].argsort()]
        #for i in np.arange(np.shape(dataNotTilted)[0]): print(dataNotTilted[i,:,0])

        print("NotTiltedShape",np.shape(dataNotTilted))
        for snapshot in range(np.shape(dataMD)[-1]):
            #datasnapshot=dataMD[:,:,snapshot]
            #dataMD[:,:,snapshot]=datasnapshot[datasnapshot[:,0].argsort()]
            dataMD[:,:,snapshot]=dataMD[:,:,snapshot][dataMD[:,0,snapshot].argsort()]
            for atom in range(np.shape(dataMD)[0]):
                if dataMD[atom,1,snapshot]==3 and dataMD[atom-1,1,snapshot]==2: #and dataMD[atom-1,0,snapshot]==47:

                    scale=1
                    # print(dataMD[atom-1,:5,snapshot])
                    # print(dataMD[atom,:5,snapshot])
                    # print(dataMD[atom+1,:5,snapshot])
                    # print(dataMD[atom+2,:5,snapshot])
                    # print("")

                    # print(dataNotTilted[atom-1,:5,0])
                    # print(dataNotTilted[atom,:5,0])
                    # print(dataNotTilted[atom+1,:5,0])
                    # print(dataNotTilted[atom+2,:5,0])
                    # print("")

                    v_Ti_Tilted=dataMD[atom-1,2:5,snapshot]
                    v_Ti_NotTilted=dataNotTilted[atom-1,2:5,0]
                    #v_Ti_NotTilted=dataMD[atom-1,2:5,0]

                    angle_x,sign_x=self.compute_angle((dataMD[atom+2,2:5,snapshot]-v_Ti_NotTilted)*scale,dataNotTilted[atom+2,2:5,0]-v_Ti_NotTilted,"x")
                    angle_y,sign_y=self.compute_angle((dataMD[atom,2:5,snapshot]-v_Ti_NotTilted)*scale,dataNotTilted[atom,2:5,0]-v_Ti_NotTilted,"y")
                    angle_z,sign_z=self.compute_angle((dataMD[atom+1,2:5,snapshot]-v_Ti_NotTilted)*scale,dataNotTilted[atom+1,2:5,0]-v_Ti_NotTilted,"z")
                    
                    #if np.abs(angle_x)>90: print(snapshot);print(dataMD[atom+2,0,snapshot]); print(dataMD[atom+2,2:5,snapshot])

                    angle_list_x=np.append(angle_list_x,angle_x*sign_x)
                    angle_list_y=np.append(angle_list_y,angle_y*sign_y)
                    angle_list_z=np.append(angle_list_z,angle_z*sign_z)


                    # angle_list_x=np.append(angle_list_x,angle_x)
                    # angle_list_y=np.append(angle_list_y,angle_y)
                    # angle_list_z=np.append(angle_list_z,angle_z)

                    displacement=(v_Ti_Tilted-v_Ti_NotTilted)
                    #displacement=np.abs(v_Ti_Tilted-v_Ti_NotTilted)
                    d_list_x=np.append(d_list_x,displacement[0])
                    d_list_y=np.append(d_list_y,displacement[1])
                    d_list_z=np.append(d_list_z,displacement[2])



                    break


        return angle_list_x,angle_list_y,angle_list_z,d_list_x,d_list_y,d_list_z




