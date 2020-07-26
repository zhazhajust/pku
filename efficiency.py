import numpy as np
import matplotlib.pyplot as plt
import sdf
import os
import scipy.signal as signal
import constant as const
import function as func
import scipy.fftpack as fftpack
import multiprocessing
from matplotlib.ticker import MultipleLocator, FuncFormatter
plt.switch_backend('agg')
limit_min=0.1e12
limit_max=10e12
locate=1800e-6
locate2=3200e-6
savefigdir=const.figdir+'Thz_'+'efficiency.png'
def draw(x):
	#p "draw",x

	savefigdir=const.figdir+str(x)+'k_bz.png'
	sdfdir=const.sdfdir +str(x).zfill(const.filenumber)+".sdf"
	if os.path.exists(sdfdir)==False:
                return [0,0]

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
		mi = 3e8/limit_min 
		ma = 3e8/limit_max
		if 2 * 3.14 / ma  > n * delta_k and  n * delta_k > 2 * 3.14 / mi:
			k_n.append(n)
	#print("n",k_n[0],k_n[-1])
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
	#print("efficiency",E_x,E_Thz,eff)

	#fig,axs=plt.subplots(2,2)
	#im=axs[0][0].pcolormesh(bz,cmap=plt.get_cmap('bwr'))	
	#im2=axs[1][0].pcolormesh(abs(k_bz),cmap=plt.get_cmap('bwr'))
	#im3=axs[0][1].pcolormesh(bz_filter.real,cmap=plt.get_cmap('bwr'))
	#im4=axs[1][1].pcolormesh(abs(k_bz2),cmap=plt.get_cmap('bwr'))
	#fig.savefig(savefigdir,dpi=200)
	#plt.close('all')
	return [E_x,E_Thz]
print(const.window_start_time)
a1 = (locate/3e8 + const.window_start_time)/const.dt_snapshot
a2 = ((locate2-800e-6)/3e8 + const.window_start_time)/const.dt_snapshot

pool = multiprocessing.Pool(processes=4)
final_energe = pool.map(draw,range(int(a1),int(a2)))
b = const.x_max/3e8/const.dt_snapshot/2
b = int(b)
print('b',b)
start=draw(b)
print('sdf1,sdf2',a1,a2)
print('Thz',limit_min,limit_max)
efficiency=np.array(final_energe/start[0])[...,1]
max_index = a1+efficiency.argmax()
max_distance = 3e8 * (max_index * const.dt_snapshot - const.window_start_time) * 1e6
print('max_index',str(max_index))
print('distance:'+str(max_distance))
print('eff:'+str(efficiency.max()))
np.savetxt(const.txtdir + 'eff.txt',efficiency)
time = np.arange(int(a1),int(a2))
locate = (time*const.dt_snapshot - const.window_start_time)*3e8*1e6
np.savetxt(const.txtdir + 'eff_locate.txt',locate)
plt.plot(locate,efficiency)
plt.savefig(savefigdir,dpi=200)
