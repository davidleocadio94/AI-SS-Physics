import ASEMDClass

STO_2x=ASEMDClass.LatticeStructureWithAtoms([1,1,1],3.8976656066670001,"../")
import ase.md
# from ase.md import MDLogger
# from ase.md import VelocityVerlet
from ase import units
dyn = ase.md.langevin.Langevin(atoms=STO_2x.OriginalCell,friction=1e-2, timestep=2*units.fs,temperature_K=300,trajectory="trj")
dyn.attach(ase.md.MDLogger(dyn, STO_2x.OriginalCell, 'md.log', header=True, stress=False,
           peratom=True, mode="a"), interval=1)
dyn.run(100000)  # take 1000 steps
