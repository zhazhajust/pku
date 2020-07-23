import numpy as np
import matplotlib.pyplot as plt
import sdf
import os
import scipy.signal as signal
import constant as const
import function as func
import scipy.fftpack as fftpack
from matplotlib.ticker import MultipleLocator, FuncFormatter
import multiprocessing
plt.switch_backend('agg')
limit_min=2e12
limit_max=5.5e12
locate=1610e-6
locate2=3000e-6
freqs = 4e12
savefigdir=const.figdir+'_'+str(freqs)+'_'+'a0.png'
def draw(x):
	c=3e8
	#w0=4e12
	me=9.10956e-31
	e=1.602176565e-19

	#p "draw",x

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

	E0_w0=bz_filter.real.max()
	print('E0:'+str(E0_w0))
	w0=freqs*2*3.1415926
	a0_w0=e*E0_w0/(me*c*w0)
	return a0_w0
#print const.window_start_time
a1 = (locate/3e8 + const.window_start_time)/const.dt_snapshot
a2 = ((locate2-800e-6)/3e8 + const.window_start_time)/const.dt_snapshot
a0 = []

index_n = 0
max_a0 = 0

pool = multiprocessing.Pool(processes=4)
#for i in range(start,stop+step,step):
#       results.append(pool.apply_async(extract, (i, ))) 
a0 = pool.map(draw,range(int(a1),int(a2)))

#for res in results:
#	a0.append(res.get())
'''
for i in range(int(a1),int(a2)):
	#print "i",i
	#eff=draw(i)
	a0_i = draw(i)
	#print eff
	a0.append(a0_i)
	if a0_i >= max_a0:
		max_a0=a0_i
		index_n = i
		#print 'eff,i',final_energe,i
'''
#print(draw(141))
#print 'sdf1,sdf2',a1,a2
#print 'w0',limit_min,limit_max
#print 'index_n',index_n
#print 'index_x',3e8 * (index_n * const.dt_snapshot - const.window_start_time) * 1e6
#print 'efficiency',max_energe/start[0]
save_a0=np.array(a0)


np.savetxt(const.txtdir + 'a0.txt',save_a0)

index_d=((a1+save_a0.argmax())*const.dt_snapshot - const.window_start_time)*3e8*1e6
print('max_a0  index   max:',index_d,save_a0.max())
time = np.arange(int(a1),int(a2))
locate = (time*const.dt_snapshot - const.window_start_time)*3e8*1e6

np.savetxt(const.txtdir + 'a0_distance.txt',locate)

#locate = (time*const.dt_snapshot - const.window_start_time)*3e8*1e6
plt.plot(locate,a0)
plt.savefig(savefigdir,dpi=200)



