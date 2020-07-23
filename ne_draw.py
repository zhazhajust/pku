import numpy as np
import matplotlib.pyplot as plt
import sdf
import os
import scipy.signal as signal
import constant as const
import imageio
import function as func
import scipy.fftpack as fftpack
from matplotlib.ticker import MultipleLocator, FuncFormatter
from multiprocessing.dummy import Pool as ThreadPool
#import wigner
plt.switch_backend('agg')
###
####
interval=50
image_list=[]
png_savedir = "./gif/png/"

def x_formatter(x, pos):
        a=(const.delta_x*x + const.c*T)   *1e6 
        return  "%d"%int(a)

def draw(x):
	#print "draw",x
	savefigdir=const.figdir+str(x)
	sdfdir=const.sdfdir +str(x).zfill(const.filenumber)+".sdf"
	data=sdf.read(sdfdir,dict=True)
	#data=sdf.read(const.sdfdir+str(x).zfill(const.filenumber)+".sdf",dict=True)
	Bz=data['Electric Field/Ey']
	global T
	time=data['Header']['time']
	if time-const.window_start_time<0:
		T=0	
	else:
		T=time-const.window_start_time	
	density=data['Derived/Number_Density/electron1'].data
	density=density.T
	ne=data['Derived/Number_Density/electron1'].data
	ne_y0=ne[...,int(const.Ny/2)]
	ne_y0=ne_y0[1900:2300]
	fig,axs=plt.subplots()
	#plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

	#fig.colorbar(im,ax=axs[0][0])
	#Xf_list=wigner.wigner(x)

	#im3=axs[1][0].pcolormesh(Xf_list,cmap=plt.get_cmap('BuPu'),rasterized=True)
	#im=axs[1][0].imshow(np.abs(zxx),extent=[],cmap=plt.get_cmap('BuPu'))
	line3=axs.plot(ne_y0,'g')
	fig.savefig(savefigdir+'ne555',dpi=200)
draw(616)
