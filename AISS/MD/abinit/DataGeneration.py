import MDDataGeneratorClass as MDDataGenerator
import numpy as np

a0=3.8976656066670001
#delta_a_list=a0*np.append(np.linspace(-.2,.2,20,endpoint=True),0)
delta_a_list=np.array([0])
for delta_a in delta_a_list:
#for delta_a in [0.0]:
	md_data=MDDataGenerator.MD_Data(SuperCellSize=[2,3,3],LatticeParameter=a0+delta_a\
,Temperatures=[1500,3000])
	md_data.GenerateFiles()
"""
md_data=MDDataGenerator.MD_Data(SuperCellSize=[2,2,2],3.8976656066670001+.5,[300,600,900,1200,1500,3000])
md_data.GenerateFiles()

md_data=MDDataGenerator.MD_Data(SuperCellSize=[2,2,2],3.8976656066670001+1,[300,600,900,1200,1500,3000])
md_data.GenerateFiles()

md_data=MDDataGenerator.MD_Data(SuperCellSize=[2,2,2],3.8976656066670001+1.5,[300,600,900,1200,1500,3000])
md_data.GenerateFiles()

md_data=MDDataGenerator.MD_Data(SuperCellSize=[2,2,2],3.8976656066670001+2.0,[300,600,900,1200,1500,3000])
md_data.GenerateFiles()"""
