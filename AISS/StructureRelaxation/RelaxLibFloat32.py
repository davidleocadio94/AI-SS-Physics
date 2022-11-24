from ase import Atoms
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.metrics import r2_score
import os
import re
from nequip.ase import nequip_calculator
from nequip.utils import Config
from ase.io import read
from ase.optimize import BFGS
from ase.build import make_supercell
import torch
import copy
import ase.md
from ase import units



#torch.set_default_dtype(torch.float64)


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
        # numatoms = File.readline()
        # h3 = File.readline()
        # xlen = File.readline()
        # ylen = File.readline()
        # zlen = File.readline()
        # h4 = File.readline()

        for a in range(NumAtoms):
            #Read string -> strip '\n' char -> split into new list
            line = File.readline().strip('\n').split()
            data[a,:,t] = line[1:4]

        t += 1

    File.close()
    return data



def z_to_number(z):
    if z=="Sr":
        return 1
    elif z=="Ti":
        return 2
    elif z=="O":
        return 3

def write_xyz(atoms_object,Type):

    Z=atoms_object.get_chemical_symbols()
    R=atoms_object.get_positions()
    f = open(f"xyz.sto.{len(Z)}.atoms.{Type}", "a")

    f.write("ITEM: TIMESTEP\n")
    f.write("0\n")
    f.write("ITEM: NUMBER OF ATOMS\n")
    f.write(f"{len(Z)}\n")
    f.write("ITEM: BOX BOUNDS pp pp pp\n")
    f.write(f"{0} {atoms_object.get_cell()[0,0]}\n")
    f.write(f"{0} {atoms_object.get_cell()[0,0]}\n")
    f.write(f"{0} {atoms_object.get_cell()[0,0]}\n")
    f.write("ITEM: ATOMS id type x y z\n")
    for i in range(1,len(Z)+1):
        z=Z[i-1]
        r=R[i-1]
        s = f"{i} {z_to_number(z)} {r[0]} {r[1]} {r[2]}\n"
        f.write(s)
    f.close()




def write_geometry(atoms_object,Type):
    L=atoms_object.get_cell()[0,0]
    Z=atoms_object.get_chemical_symbols()
    R=atoms_object.get_positions()
    f = open(f"geometry.{len(Z)}.in.{Type}", "a")
    f.write(f"lattice_vector {L} {0} {0}\n")
    f.write(f"lattice_vector {0} {L} {0}\n")
    f.write(f"lattice_vector {0} {0} {L}\n")

    for i in range(len(Z)):
        z=Z[i]
        r=R[i]/L
        f.write(f"atom_frac {r[0]} {r[1]} {r[2]} {z}\n")

    f.close()


def read_dump_log_lammps(FileName,Nfield, SnapShots):
        """read lammps dumpfile with header (not saved)"""

        File = open(FileName,'r')
        data = np.ndarray((SnapShots,Nfield),dtype=float)
        
        t = 0
        line_exit=0
        while(line_exit==0):
            h=File.readline()
            #print(h)
            if "Time Step Temp PotEng TotEng Press Volume " in h:
                line_exit=1

        while (t < SnapShots):
            #Read string -> strip '\n' char -> split into new list
            line = File.readline().strip('\n').split()
            data[t,:] = line
                
            t += 1

        File.close()
        return data

def read_dump_geometry_xyz(tmin,FileName,NumAtoms,Nfield,SnapShots):
    """read lammps dumpfile with header (not saved)"""

    File = open(FileName,'r')
    data = np.ndarray((NumAtoms,Nfield),dtype=float)
    
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
            if t==tmin:
                data[a,:] = line


            
        t += 1
    data=data[data[:,0].argsort()]
    File.close()
    return data









class LatticeStructureWithAtoms():

    def __init__(self,SS_Size,a,ModelPath,ftol=.01):


        self.a=a
        self.SS_Size=SS_Size
        self.ModelPath=ModelPath
        self.ftol=ftol

        NN=nequip_calculator.NequIPCalculator
        self.NN=NN.from_deployed_model(model_path=ModelPath, \
       species_to_type_name = {
           "Sr": "Sr",
           "Ti":"Ti",
           "O":"O"
       })



        self.OriginalCell=Atoms('SrTiO3',positions=[(0, 0, 0),
                                    (0.5*a, 0.5*a, 0.5*a),
                                    (0.0*a, 0.5*a, 0.5*a),
                                    (0.5*a, 0*a, 0.5*a),
                                    (0.5*a, 0.5*a, 0*a)],\
                                        cell=[1*a,1*a,1*a],pbc=True)
        multiplier = np.identity(3) * SS_Size
        self.OriginalCell = make_supercell(self.OriginalCell, multiplier)
        self.OriginalCell.set_calculator(self.NN)

        
        self.OriginalEnergy=self.OriginalCell.get_potential_energy()
        #self.OriginalPositions=self.OriginalCell.get_positions()
        self.OriginalLatticeParameter=a

        self.NewestCell=copy.deepcopy(self.OriginalCell)
        # print("NewestCell id",id(self.NewestCell))
        # print("OriginalCell id",id(self.OriginalCell))

        self.NewestEnergy=self.NewestCell.get_potential_energy()
        self.NewestPositions=self.NewestCell.get_positions().copy()
        self.NewestLatticeParameter=a 
        self.NewestCellcell=self.NewestCell.cell.copy()

        self.NAtoms=5*self.SS_Size**3   

        self.EnergyDropList=[]
        self.EnergyList=[]
        self.LatticeParameterList=[]
        # print("NewestPositions id",id(self.NewestPositions))
        # print("NewestCellgetpositions id",id(self.NewestCell.positions))      
    def RunMD(self,dt,steps,temp):

        dyn = ase.md.langevin.Langevin(atoms=self.NewestCell,friction=1e-2, timestep=dt*units.fs,temperature_K=temp,trajectory="trj")
        dyn.attach(ase.md.MDLogger(dyn, self.NewestCell, 'md.log', header=True, stress=False, peratom=True, mode="a"), interval=1)
        dyn.run(steps)  # take 1000 steps    

        #self.OriginalEnergy=self.OriginalCell.get_potential_energy()

        #self.NewestCell=copy.deepcopy(self.OriginalCell)


        self.NewestEnergy=self.NewestCell.get_potential_energy()
        self.NewestPositions=self.NewestCell.get_positions().copy()


    def UpdateCellWithMDMin(self,PathToMDDump,PathToXYZDump,NumAtoms):

        data = read_dump_log_lammps(PathToMDDump,7,500)

        data=data[data[:,3].argsort()]
        tmin=int(data[0,1])
        data=read_dump_geometry_xyz(tmin,PathToXYZDump,NumAtoms,8,50000)[:,2:5]
        self.NewestCell.positions=data

        self.NewestEnergy=self.NewestCell.get_potential_energy()
        self.NewestPositions=self.NewestCell.get_positions().copy()
        self.NewestCellcell=self.NewestCell.cell.copy()


    def ComputeRelaxationOfAtomicPositions(self,iteration):

        self.NewestCell.positions=(self.NewestCell.positions+\
        self.NewestCell.positions*np.random.normal(0,.0005/(iteration+1),np.shape(self.NewestCell.positions)))
        opt = BFGS(self.NewestCell)
        opt.run(fmax=self.ftol)
        print("lowering newest cell with no lattice\
             relaxation ",(self.NewestCell.get_potential_energy()\
                 -self.OriginalEnergy)/(self.SS_Size**3))
        self.NewestPositions=self.NewestCell.positions.copy()
 
    def ComputeRelaxationOfLattice(self):

        print("Checking if it's relaxed ",(self.NewestCell.get_potential_energy()\
            -self.OriginalEnergy)/(self.SS_Size**3))

        a_list=self.SS_Size*(self.NewestLatticeParameter+np.linspace(-.1,.1,20))
        Energy_list=np.array([])
        for a in a_list:
            # print("a",a)
            # print("olda",self.NewestLatticeParameter)
            self.NewestCell.positions=self.NewestPositions.copy()*a/\
                (self.SS_Size*self.NewestLatticeParameter)
            self.NewestCell.cell=self.NewestCellcell.copy()*a/\
              (self.SS_Size*self.NewestLatticeParameter)
            Energy_list=np.append(Energy_list,self.NewestCell.get_potential_energy())
        
        index_min = np.argmin(Energy_list)
        amin=a_list[index_min]/self.SS_Size
        print(f"a = {amin}, Emin= {Energy_list[index_min]}")
        Emin=Energy_list[index_min]
        Edrop=(Emin-self.OriginalEnergy)/(self.SS_Size**3)
        print(f"Eorig={self.OriginalEnergy}")
        print(f"Edrop is {Edrop}")

        
        self.NewestPositions=self.NewestPositions*amin/\
                (self.NewestLatticeParameter)
        self.NewestCellcell=self.NewestCellcell*amin/\
                (self.NewestLatticeParameter)
        self.NewestLatticeParameter=amin
        self.NewestEnergy=self.NewestCell.get_potential_energy()

        self.NewestCell.positions=self.NewestPositions.copy()
        self.NewestCell.cell=self.NewestCellcell.copy()
        # print(type(self.NewestCell.positions))
        # print(id(self.NewestCell.positions))
        # self.NewestCell.positions=self.NewestPositions*amin/\
        #         (self.NewestLatticeParameter)
        # self.NewestCell.cell=self.NewestCellcell*amin/\
        #         (self.NewestLatticeParameter)
        self.NewestEnergy=self.NewestCell.get_potential_energy()       


        self.EnergyDropList.append(Edrop)
        self.EnergyList.append(self.NewestEnergy)
        self.LatticeParameterList.append(amin)
        
