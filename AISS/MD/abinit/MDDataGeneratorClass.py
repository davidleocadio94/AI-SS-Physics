import os
import os.path
class MD_Data():
    def __init__(self,SuperCellSize,LatticeParameter,Temperatures):
        self.SuperCellSize=SuperCellSize
        self.LatticeParameter=LatticeParameter
        self.Temperatures=Temperatures
        self.GrandFatherPath=os.getcwd()
        self.ParentPath="MD_SS_"+str(self.SuperCellSize).replace(" ","").replace("[","").replace("]","").replace(",","x")+\
"_LattParam_"+str(self.LatticeParameter)
        self.ParentPath=os.path.join(os.getcwd(),self.ParentPath)
        print(self.ParentPath)

    def GenerateFiles(self):
        if os.path.isdir(self.ParentPath):
            print("File already exists, no new directory created")
        else:
            os.makedirs(self.ParentPath)
            print("File created")
        os.chdir(self.ParentPath)

        PrimitiveGeometryText="lattice_vector "+ str(f"{self.LatticeParameter:.8f}") + " 0.0000000000000000 0.0000000000000000 \n \
lattice_vector 0.0000000000000000 "+str(f"{self.LatticeParameter:.8f}") +" 0.0000000000000000 \n \
lattice_vector 0.0000000000000000 0.0000000000000000 " + str(f"{self.LatticeParameter:.8f}")+ " \n \
atom_frac -0.0000002400000000 -0.0000005700000000 -0.0000000300000000 Sr \n \
atom_frac 0.4999997600000000 0.4999994300000000 0.4999999700000000 Ti \n \
atom_frac -0.0000002400000000 0.4999994300000000 0.4999999700000000 O \n \
atom_frac 0.4999997600000000 -0.0000005700000000 0.4999999700000000 O \n \
atom_frac 0.4999997600000000 0.4999994300000000 -0.0000000300000000 O"

        GeometryPrimitiveFileName="geometry.in.primitive"

        if os.path.isfile(GeometryPrimitiveFileName):
            print("Primitive geometry file exists, no new creation")
        else:
            with open(GeometryPrimitiveFileName,'w') as f: f.write(PrimitiveGeometryText)
            print("File didnt exist so primitive geometry file was created")
	
        for Temperature in self.Temperatures:
            os.chdir(self.ParentPath)
            PathForTemperature="MD_"+str(Temperature)
            if os.path.isdir(PathForTemperature):
                print("There is data at this temperature and lattice param already")
            else:
                print("There is no data at this temperature, proceeding to generating data")
                os.makedirs(PathForTemperature)
                os.chdir(PathForTemperature)
                os.system("ln -s ../geometry.in.primitive")
                os.system("vibes utils make-supercell geometry.in.primitive -dd "+\
                    str(self.SuperCellSize[0])+" "+str(self.SuperCellSize[1])+ " "+str(self.SuperCellSize[2]))
                os.system("mv geometry.in.primitive.supercell* geometry.in.supercell")
                os.system("vibes utils create-samples geometry.in.supercell -T "+str(Temperature))
                os.system("mv geometry.in.supercell.* geometry.in")

                MdInFile="""[machine]
basissetloc:                   /u/vdavi/FHIaims/species_defaults
aims_command:                  /u/vdavi/FHIaims/run_aims.sh

[files]
geometry:                      geometry.in
primitive:                     geometry.in.primitive
supercell:                     geometry.in.supercell

[calculator]
name:                          aims
socketio:                      True

[calculator.parameters]
xc:                            pbesol
charge_mix_param:              0.1
compute_forces:                True
sc_accuracy_rho:               1e-06
relativistic:                  atomic_zora scalar
output_level:                  MD_light
use_pimd_wrapper:              ('localhost', 10011)

[calculator.kpoints]
density:                       5

[calculator.basissets]
Sr:                            light
Ti:                            light
O:                             light

[md]
driver:                        Langevin
timestep:                      5
maxsteps:                      10000
#[restart]
#command = sbatch submit.sh
compute_stresses:              False
workdir:                       md

[md.kwargs]
temperature:                   """+str(Temperature)+"""
friction:                      0.02
logfile:                       log.md"""
        
            with open('md.in','w') as f: f.write(MdInFile)

            myprog="""#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.out.%j
#SBATCH -e ./tjob.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J """+str(self.SuperCellSize).replace(" ","").replace("[","").replace("]","").replace(",","x")+str(Temperature)+"K"+"""
#
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=40
#
#SBATCH --mail-type=none
#SBATCH --mail-user=userid@example.mpg.de
#
# Wall clock limit:
#SBATCH --time=23:00:00

#module purge
#module load intel/18.0.5 impi/2018.4 mkl/2018.4 anaconda/3/2020.02
#rm tjob*
#rm -rf md
vibes run md md.in &> log.md"""

            with open('myprog','w') as f: f.write(myprog)     
            os.system("sbatch ./myprog > prog.out")
        os.chdir(self.GrandFatherPath)
