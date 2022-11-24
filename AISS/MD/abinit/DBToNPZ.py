    


from ase.io import read
import numpy as np
import os

ParentPath=os.getcwd()

path_list=["MD_SS_3x3x3_LattParam_4.897665606667","MD_SS_3x3x3_LattParam_5.897665606667","MD_SS_3x3x3_LattParam_4.397665606667","MD_SS_3x3x3_LattParam_5.397665606667","MD_SS_3x3x3_LattParam_8.897665606667001"]



for path_data in path_list:
 for temperature in ["MD_1500","MD_300"]:
        os.chdir(ParentPath)
        os.chdir(path_data+"/"+temperature)
        present=os.getcwd(); print(present) 
        os.system("vibes utils trajectory 2db md/trajectory.son")
        atoms = read('trajectory.db',index=':')
        
        
        # parse properties as list of dictionaries
        property_list = []
        E=[]
        F=[]
        R=[]
        z=[]
        CELL=[]
        PBC=[]
        for at in atoms:
        
            E.append([at.get_potential_energy()])
            F.append(at.get_forces())
            R.append(at.get_positions())
            z.append(at.numbers)
            CELL.append(at.get_cell())
            PBC.append(at.get_pbc())
        
        E=np.array(E)
        F=np.array(F)
        R=np.array(R)
        z=np.array(z)
        CELL=np.array(CELL)
        PBC=np.array(PBC)
        #randomly_permuted_idx=np.random.permutation(np.arange(len(E)))
        randomly_permuted_idx=np.arange(len(E))
        size_test=1
        #randomly_permuted_idx=np.arange(len(atoms))
        idxsTrain=randomly_permuted_idx[:-size_test]
        idxsTrain=idxsTrain.astype(int)
        idxsTest=randomly_permuted_idx[-size_test:]
        idxsTest=idxsTest.astype(int)
        #idxMINS = np.argpartition(E.flatten(), 5)
#        print(idxMINS)
        #idxsTest=np.append(idxsTest,idxMINS[:5])
        len(idxsTest)+len(idxsTrain)
        
        
        np.savez('trajectory',E=E[idxsTrain],F=F[idxsTrain],R=R[idxsTrain],z=z[0],CELL=\
                CELL[0],PBC=PBC[0])
        np.savez('trajectoryTest',E=E[idxsTest],F=F[idxsTest],R=R[idxsTest],z=z[0],CELL=\
                CELL[0],PBC=PBC[0])
        print("Train",np.shape(idxsTrain))
        print("Test",np.shape(idxsTest))
