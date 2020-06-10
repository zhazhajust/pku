import numpy as np
import matplotlib.pyplot as plt
import sdf
import os
import scipy.signal as signal
import constant as const
import function as func
import scipy.fftpack as fftpack
from matplotlib.ticker import MultipleLocator, FuncFormatter
plt.switch_backend('agg')
import constant as const
delta_k=3.14/const.delta_x/(const.Nx/2)
k_n=[]
for n in range(0,const.Nx):
	if 2 * 3.14 / 5e-6  > n * delta_k and  n * delta_k > 2 * 3.14 / 150e-6:
		k_n.append(n)
print(k_n)
def draw(x):
        #p "draw",x
        savefigdir=const.figdir+str(x)+'k_bz.png'
        sdfdir=const.sdfdir +str(x).zfill(const.filenumber)+".sdf"
        data=sdf.read(sdfdir,dict=True)
        Bz=data['Electric Field/Ey']
        time=data['Header']['time']
        bz=Bz.data
        density=data['Derived/Number_Density/electron1'].data
        bz=bz.T
        k_bz=np.fft.fft(bz)
	bz_f=np.fft.ifft(k_bz)
	print(bz,bz_f.real)
	print(k_bz)
draw(4000)
