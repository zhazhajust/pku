import numpy as np
import constant as const
import sdf
savefigdir=const.figdir+str(2000)+'k_bz.png'
sdfdir=const.sdfdir +str(2000).zfill(const.filenumber)+".sdf"
data=sdf.read(sdfdir,dict=True)
Bz=data['Electric Field/Ey']
np.savetxt('aaa.txt',Bz)
