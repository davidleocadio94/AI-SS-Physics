import RelaxLib

SS_List=[6,8]
iterations=15


for SS in SS_List:
      
     STO=RelaxLib.LatticeStructureWithAtoms(SS,3.8976656066670001,".",ftol=0.01)
     #STO.UpdateCellWithMDMin("/u/vdavi/Documents/135_atoms/BIGDATASETS/IterativeRelaxationEdrop/OOP/15KDataFromRaven/log.lammps",
#"/u/vdavi/Documents/135_atoms/BIGDATASETS/IterativeRelaxationEdrop/OOP/15KDataFromRaven/dump.velocity",int(5*SS**3))
     
     STO.RunMD(dt=2,steps=200,temp=15)

     for i in range(iterations):
         STO.ComputeRelaxationOfAtomicPositions(iteration=i)
         STO.ComputeRelaxationOfLattice()
     print("SuperCell Size",STO.SS_Size)
     print("aoriginal ",STO.OriginalLatticeParameter)
     print("anew",STO.NewestLatticeParameter)
     print("Eold",STO.OriginalEnergy)
     print("Enew",STO.NewestEnergy)
     print("alist",STO.LatticeParameterList)
     print("EDropList",STO.EnergyDropList)
     print("EList",STO.EnergyList)
     
     RelaxLib.write_geometry(STO.OriginalCell,"Cubic")
     RelaxLib.write_geometry(STO.NewestCell,"Relaxed")
     
     
     
     RelaxLib.write_xyz(STO.OriginalCell,"Cubic")
     RelaxLib.write_xyz(STO.NewestCell,"Relaxed")
