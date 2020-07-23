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
limit_min=0.1e12
limit_max=10e12
locate=8500e-6
locate2=10000e-6
def draw(x):
        #p "draw",x
        savefigdir=const.figdir+str(x)+'k_bz.png'
        sdfdir=const.sdfdir +str(x).zfill(const.filenumber)+".sdf"
        data=sdf.read(sdfdir,dict=True)
        Bz=data['Electric Field/Ey']
	time=data['Header']['time']
	bz=Bz.data
	bz=bz.T
	k_bz=np.fft.fft(bz)	
	delta_k=3.14/const.delta_x/(const.Nx/2)
	k_bz2=k_bz*1
	k_n=[]
	for n in range(0,const.Nx):
		mi = 3e8/limit_min 
		ma = 3e8/limit_max
		if 2 * 3.14 / ma  > n * delta_k and  n * delta_k > 2 * 3.14 / mi:
			k_n.append(n)
	k_bz2[...,0:k_n[0]]=0    #k_bz.argmin()
	k_bz2[...,k_n[-1]:-k_n[-1]]=0  #k_bz.argmin()
	k_bz2[...,-k_n[0]:]=0    #k_bz.argmin()
	bz_filter=np.fft.ifft(k_bz2)
	E_x=np.sum(np.sum(np.square(bz)))
	E_Thz=np.sum(np.sum(np.square(bz_filter.real)))
	eff=E_Thz/E_x
	return [E_x,E_Thz]
print const.window_start_time
middle = (locate/3e8 + const.window_start_time)/const.dt_snapshot
middle = int(middle)
a1 = (locate/3e8 + const.window_start_time)/const.dt_snapshot
a2 = ((locate2-800e-6)/3e8 + const.window_start_time)/const.dt_snapshot
d_n = 800e-6/3e8/const.dt_snapshot
d_n = int(d_n)
#print 'middle',middle
#print "d_n",d_n
final_energe = []
index_n = 0
max_energe = 0
#for i in range(middle - int(d_n),middle + d_n/8,1):
print 'a1,a2',a1,a2
for i in range(int(a1),int(a2)):
	print "i",i
	eff=draw(i)
	print eff
	final_energe.append(eff[1])
	if eff[1] >= max_energe:
                max_energe=eff[1]
                index_n = i
                #print 'eff,i',final_energe,i
#final=draw(a)
b = const.x_max/3e8/const.dt_snapshot/2
b = int(b)
print 'b',b
start=draw(b)
print 'sdf1,sdf2',a1,a2
print 'Thz',limit_min,limit_max
print 'index_n',index_n
print 'index_x',3e8 * (index_n * const.dt_snapshot - const.window_start_time) * 1e6
print 'efficiency',max_energe/start[0]
efficiency=np.array(final_energe/start[0])
np.savetxt(const.txtdir + 'eff.txt',efficiency)
