
import MDDistLib
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kurtosis
from scipy.stats import skew
import pandas as pd
import seaborn as sns
import sys
from scipy.stats import gaussian_kde

natoms=int(float(sys.argv[1]))
Temp=sys.argv[2]

MD_4x4x4=MDDistLib.MDFile(PathMD=str(Temp)+"/dump.velocity",\
    PathNotTilted="/home/villarreal/Documents/NN-MD/LAMMPSAndGeometryAnalysis/OOP/xyz.sto.320.atoms.Cubic",\
        NumAtoms=natoms,TotalSnapShotsForAvg=50000,TotalMDLength=50000)


angle_list_x,angle_list_y,angle_list_z,d_list_x,d_list_y,d_list_z=MD_4x4x4.GetAngleDistMD()

angle_list_x = [i for i in angle_list_x if np.abs(i)<90]
angle_list_y = [i for i in angle_list_y if np.abs(i)<90]
angle_list_z = [i for i in angle_list_z if np.abs(i)<90]


length=np.min([len(angle_list_x),len(angle_list_y),len(angle_list_z)])

angle_list_x=angle_list_x[:length];angle_list_y=angle_list_y[:length]; angle_list_z=angle_list_z[:length]

print("Temperature statistics for "+Temp+"K")

scipy_kernel=gaussian_kde(np.abs(angle_list_x))
u=np.linspace(np.min(angle_list_x),np.max(angle_list_x),1000)
v=scipy_kernel.evaluate(u)
IdxMax=np.argmax(v)
print("ThetaX ",u[IdxMax])

scipy_kernel=gaussian_kde(np.abs(angle_list_y))
u=np.linspace(np.min(angle_list_y),np.max(angle_list_y),1000)
v=scipy_kernel.evaluate(u)
IdxMax=np.argmax(v)
print("ThetaY ",u[IdxMax])

scipy_kernel=gaussian_kde(np.abs(angle_list_z))
u=np.linspace(np.min(angle_list_z),np.max(angle_list_z),1000)
v=scipy_kernel.evaluate(u)
IdxMax=np.argmax(v)
print("ThetaZ ",u[IdxMax])

#print("Temperature statistics for "+Temp+"K")
#
#counts,bin_edges=np.histogram(np.abs(angle_list_x))
#print("ThetaXMax = ",bin_edges[np.argmax(counts)])
#
#
#counts,bin_edges=np.histogram(np.abs(angle_list_y))
#print("ThetaYMax = ",bin_edges[np.argmax(counts)])
#
#counts,bin_edges=np.histogram(np.abs(angle_list_z))
#print("ThetaZMax = ",bin_edges[np.argmax(counts)])

plt.figure()

df = pd.DataFrame({'thetax':angle_list_x, 'thetay':angle_list_y, 'thetaz':angle_list_z})


df.columns.name="AngleDirection"
df=df.stack()
df.name="AngleValues (deg)"
df=df.reset_index()


g=sns.histplot(
    df, x="AngleValues (deg)", hue="AngleDirection", element="step",
    stat="density", common_norm=False,multiple="stack"
)
plt.legend(labels=["ThetaX","ThetaY","ThetaZ"],bbox_to_anchor=(1.04,1), loc="upper left")
plt.savefig(str(Temp)+"AngleDistSNS.png",bbox_inches='tight')



plt.figure()


df = pd.DataFrame({'DeltaX':d_list_x, 'DeltaY':d_list_y, 'DeltaZ':d_list_z})


df.columns.name="TiDispDirection"
df=df.stack()
df.name="DeltaValues (Ang)"
df=df.reset_index()



g=sns.histplot(
    df, x="DeltaValues (Ang)", hue="TiDispDirection", element="step",
    stat="density", common_norm=False,multiple="stack"
)
plt.legend(labels=["DeltaX","DeltaY","DeltaZ"],bbox_to_anchor=(1.04,1), loc="upper left")
plt.savefig(str(Temp)+"TiDispSNS.png",bbox_inches='tight')






###plt.figure()
##bins = np.linspace(-30, 30, 100)
##_=plt.hist(angle_list_x,bins,alpha=0.25,density=True,label=f"Theta_x\n avg={np.mean(angle_list_x):.2f} std = {np.std(angle_list_x):.2f} skw = {skew(angle_list_x):.2f} krt = {kurtosis(angle_list_x):.2f} ")
##_=plt.hist(angle_list_y,bins,alpha=0.25,density=True,label=f"Theta_y\n avg={np.mean(angle_list_y):.2f} std = {np.std(angle_list_y):.2f} skw = {skew(angle_list_y):.2f} krt = {kurtosis(angle_list_y):.2f}")
##_=plt.hist(angle_list_z,bins,alpha=0.25,density=True,label=f"Theta_z\n avg={np.mean(angle_list_z):.2f} std = {np.std(angle_list_z):.2f} skw = {skew(angle_list_z):.2f} krt = {kurtosis(angle_list_z):.2f}")
###plt.hist(theta_list_lowT*180/np.pi,bins,alpha=0.5,label='185.0K')
##plt.ylabel("Octahedral tilting population")
##plt.xlabel("Angle (degrees)")
### plt.hist(x, bins, alpha=0.5, label='x')
### plt.hist(y, bins, alpha=0.5, label='y')
##plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
##plt.savefig(str(Temp)+"xyzAngleDist.png",bbox_inches='tight')
##
##
##
##
##plt.figure()
##bins = np.linspace(-1, 1, 100)
##_=plt.hist(d_list_x,bins,alpha=0.25,density=True,label=f"d_x\n avg={np.mean(d_list_x):.2f} std = {np.std(d_list_x):.2f} skw = {skew(d_list_x):.2f} krt = {kurtosis(d_list_x):.2f} ")
##_=plt.hist(d_list_y,bins,alpha=0.25,density=True,label=f"d_y\n avg={np.mean(d_list_y):.2f} std = {np.std(d_list_y):.2f} skw = {skew(d_list_y):.2f} krt = {kurtosis(d_list_y):.2f}")
##_=plt.hist(d_list_z,bins,alpha=0.25,density=True,label=f"d_z\n avg={np.mean(d_list_z):.2f} std = {np.std(d_list_z):.2f} skw = {skew(d_list_z):.2f} krt = {kurtosis(d_list_z):.2f}")
###plt.hist(theta_list_lowT*180/np.pi,bins,alpha=0.5,label='185.0K')
##plt.ylabel("Ti disp population")
##plt.xlabel("displacement (Angs)")
### plt.hist(x, bins, alpha=0.5, label='x')
### plt.hist(y, bins, alpha=0.5, label='y')
##plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
##plt.savefig(str(Temp)+"xyzDispDist.png",bbox_inches='tight')
##
##
##df=pd.DataFrame({"ThetaX":angle_list_x[:10000],"ThetaY":angle_list_y[:10000]})
##sns.set(font_scale=1.2)
##g=sns.jointplot(data=df,x="ThetaX",y="ThetaY",height=6)
##g.savefig(str(Temp)+"ThetaXY.png")
##
##
##df=pd.DataFrame({"ThetaX":angle_list_x[:10000],"ThetaZ":angle_list_z[:10000]})
##sns.set(font_scale=1.2)
##g=sns.jointplot(data=df,x="ThetaX",y="ThetaZ",height=6)
##g.savefig(str(Temp)+"ThetaXZ.png")
##
##df=pd.DataFrame({"ThetaY":angle_list_y[:10000],"ThetaZ":angle_list_z[:10000]})
##sns.set(font_scale=1.2)
##g=sns.jointplot(data=df,x="ThetaY",y="ThetaZ",height=6)
##g.savefig(str(Temp)+"ThetaYZ.png")
##
###plt.style.use('seaborn-deep')
##
### #plt.figure()
### fig, axes = plt.subplots(1, 3)
### bins = np.linspace(-30, 30, 100)
### axes[0].hist(angle_list_x,bins,alpha=0.5,density=True,label=f"Theta_x, avg={np.mean(angle_list_x)}")
### axes[1].hist(angle_list_y,bins,alpha=0.5,density=True,label=f"Theta_y, avg={np.mean(angle_list_y)}")
### axes[2].hist(angle_list_z,bins,alpha=0.5,density=True,label=f"Theta_z, avg={np.mean(angle_list_z)}")
### #plt.hist(theta_list_lowT*180/np.pi,bins,alpha=0.5,label='185.0K')
### plt.ylabel("Octahedral tilting population")
### plt.xlabel("Angle (degrees)")
### # plt.hist(x, bins, alpha=0.5, label='x')
### # plt.hist(y, bins, alpha=0.5, label='y')
### plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
### #plt.savefig("1585vs185AngleDist.png")
##
##
##
