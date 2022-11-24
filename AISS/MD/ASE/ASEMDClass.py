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

torch.set_default_dtype(torch.float64)




class LatticeStructureWithAtoms():

    def __init__(self,SS_Size,a,ModelPath):


        self.a=a
        self.SS_Size=SS_Size
        self.ModelPath=ModelPath


        NN=nequip_calculator.NequIPCalculator
        self.NN=NN.from_deployed_model(model_path=ModelPath+"/deployeddouble.pth", \
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
        #multiplier = np.identity(3) * SS_Size
        multiplier=np.diag(self.SS_Size)
        self.OriginalCell = make_supercell(self.OriginalCell, multiplier)
        self.OriginalCell.set_calculator(self.NN)


        self.OriginalEnergy=self.OriginalCell.get_potential_energy()
        #self.OriginalPositions=self.OriginalCell.get_positions()
        self.OriginalLatticeParameter=a



        self.NAtoms=5*np.linalg.det(np.diag(self.SS_Size))
