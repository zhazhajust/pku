import numpy as np
import matplotlib.pyplot as plt
import sdf
plt.switch_backend('agg')
x=18000
sdfdir="../Data/acceleration/"
savefigdir="./acceleration/fig/"+str(x)
data=sdf.read(sdfdir+str(x)+".sdf",dict=True)
data
Bz=data['Magnetic Field/Bz']
bz=Bz.data
density=data['Derived/Number_Density/electron'].data
bz=bz.T
dengsity=density.T
fig,ax=plt.subplots()
fig2,ax2=plt.subplots()
im=ax.pcolormesh(bz,cmap=plt.get_cmap('rainbow'))
im2=ax2.pcolormesh(density,cmap=plt.get_cmap('gray'))
fig.savefig(savefigdir+"bz.png",dpi=200)

fig2.savefig(savefigdir+"density.png",dpi=200)
