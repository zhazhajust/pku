import numpy as np
import matplotlib.pyplot as plt
import sdf
import matplotlib.pyplot as pl
from matplotlib.ticker import MultipleLocator, FuncFormatter
plt.switch_backend('agg')
micron = 1e-6
c=3e8

locate=500 * micron
locate2=16500 * micron
dt=3e-15
t = int((locate/c)/dt)
t2 = int((locate/c)/dt)


sdfdir="../Data/density_120%/"+str(t).zfill(5)+".sdf"
sdfdir2="../Data/density_120%/"+str(t2).zfill(5)+".sdf"
savedir="./fig/density1.2/"+str(locate/micron)+".png"
data=sdf.read(sdfdir,dict=True)
Bz=data["Magnetic Field/Bz"].data


fig,ax=plt.subplots()
bz=Bz.T
im=ax.pcolormesh(bz,cmap=plt.get_cmap("rainbow"))

fig.savefig(savedir,dpi=200)
