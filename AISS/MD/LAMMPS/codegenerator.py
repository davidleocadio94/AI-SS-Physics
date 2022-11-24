import numpy as np
import os
import subprocess
#Ts=[100,  1200,  150,  1500,  300,  50,  600,  900]



#Ts=np.round(np.logspace(2.9,3.1,num=4))
#Ts=#np.round(np.logspace(1.699,3.1,num=5))
#Ts=np.round(np.logspace(1.8,3.2,num=4))
Ts=[300]
ParentPath=os.getcwd()

for T in Ts:
    os.mkdir(str(T))

    submissionscript="""#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./job.out.%j
#SBATCH -e ./job.err.%j
# Initial working directory:
#SBATCH -D ./
# Job name
#SBATCH -J +"""+str(T)+"""
#
#SBATCH --ntasks=1
#SBATCH --constraint="gpu"
#
# --- default case: use a single GPU on a shared node ---
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=18
#SBATCH --mem=125000
#
# --- uncomment to use 2 GPUs on a shared node ---
# #SBATCH --gres=gpu:a100:2
# #SBATCH --cpus-per-task=36
# #SBATCH --mem=250000
#
# --- uncomment to use 4 GPUs on a full node ---
# #SBATCH --gres=gpu:a100:4
# #SBATCH --cpus-per-task=72
# #SBATCH --mem=500000
#
#SBATCH --mail-type=none
#SBATCH --mail-user=userid@example.mpg.de
#SBATCH --time=22:30:00


module purge
module load intel/21.3.0
module load impi/2021.3
module load mkl/2021.3
module load cuda/11.4
module load cudnn/8.2.4
module load git/2.31
module load anaconda/3/2021.05
module load gcc/11
module load cmake


export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

srun ~/lammps/build/lmp -in in.sto"""










    text="""boundary p p p
units metal
atom_style atomic
newton off
#read_data structure.data
variable a equal  3.89766561



lattice custom    1     &
        a1      $a      0.0     0.0     &
        a2      0.0     $a      0.0     &
        a3      0.0     0.0     $a      &
        basis   0.0     0.0     0.0     &
        basis   0.5     0.5     0.5     &
        basis   0.0     0.5     0.5     &
        basis   0.5     0.0     0.5     &
        basis   0.5     0.5     0.0


#region myreg block 0 $a 0 $a 0 $a
region myreg block 0 4 0 4 0 4
create_box 3 myreg
create_atoms 3 box &
basis 1 1 &
basis 2 2 &
basis 3 3 &
basis 4 3 &
basis 5 3


pair_style nequip
pair_coeff * * ../deployed.pth Sr Ti O


mass 1 87.62
mass 2 47.88
mass 3 16.00

neighbor        2.0 bin
neigh_modify    delay 0 every 10 check yes


compute myRDF all rdf 50 1 1 2 2 3 3 1 2 1 3 2 3
fix 1 all ave/time 100 1 100 c_myRDF[*] file tmp.rdf mode vector





thermo          100
thermo_style    custom time step temp pe etotal press vol
thermo_modify   norm no flush yes

dump 1 all custom 1 dump.velocity id type x y z vx vy vz



velocity all create """+str(T)+""" 12345
fix             ensemble_set all nvt temp """+str(T)+""" """+str(T)+""" $(100*dt)


#write_data      waterinfo.data
#write_restart   waterinfo.restart

timestep        0.001
run             50000"""
    with open(str(T)+'/in.sto','w') as f: f.write(text) 
    with open(str(T)+'/myprog_raven','w') as f: f.write(submissionscript) 
    os.chdir(str(T))
    subprocess.run(["sbatch","./myprog_raven",">","prog.out"])
    os.chdir(ParentPath)
