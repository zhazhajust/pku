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
def draw(x):
        #p "draw",x
        savefigdir=const.figdir+str(x)+'k_bz.png'
        sdfdir=const.sdfdir +str(x).zfill(const.filenumber)+".sdf"
        data=sdf.read(sdfdir,dict=True)
        Bz=data['Electric Field/Ey']
        time=data['Header']['time']
        bz=Bz.data
        #density=data['Derived/Number_Density/electron1'].data
        bz=bz.T
        k_bz=np.fft.fft(bz)
        #from scipy import signal
        #k_min=2*const.delta_x/150e-6
        #k_max=2*const.delta_x/75e-6
        #b, a = signal.butter(8, [k_min,k_max], 'bandpass')
        #filtedData = signal.filtfilt(b, a, bz)
        # frequency
        delta_k=3.14/const.delta_x/(const.Nx/2)
        k_bz2=k_bz*1
        k_n=[]
        for n in range(0,const.Nx):
                mi = 3e8/2e12
                ma = 3e8/5e12
                if 2 * 3.14 / ma  > n * delta_k and  n * delta_k > 2 * 3.14 / mi:
                        k_n.append(n)
        print("n",k_n[0],k_n[-1])
        k_bz2[...,0:k_n[0]]=0    #k_bz.argmin()
        k_bz2[...,k_n[-1]:-k_n[-1]]=0  #k_bz.argmin()
        k_bz2[...,-k_n[0]:]=0    #k_bz.argmin()
        #k_bz1=k_bz*1
        #k_bz1[...,200:-200]=0
        #bz_filter1=np.fft.ifft(k_bz1)
        bz_filter=np.fft.ifft(k_bz2)
        E_x=np.sum(np.sum(np.square(bz)))
        E_Thz=np.sum(np.sum(np.square(bz_filter.real)))
        eff=E_Thz/E_x
        print("efficiency",E_x,E_Thz,eff)
        fig,axs=plt.subplots(2,2)
        im=axs[0][0].pcolormesh(bz,cmap=plt.get_cmap('bwr'))   
        im2=axs[1][0].pcolormesh(abs(k_bz),cmap=plt.get_cmap('bwr'))
        im3=axs[0][1].pcolormesh(bz_filter.real,cmap=plt.get_cmap('bwr'))
        im4=axs[1][1].pcolormesh(abs(k_bz2),cmap=plt.get_cmap('bwr'))
        fig.savefig(savefigdir,dpi=200)
	plt.clf()
        plt.close('all')
middle = (868e-6/3e8 + const.window_start_time)/const.dt_snapshot
middle = int(middle)
draw(middle)