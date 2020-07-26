import numpy as np
import constant as const
import os
import multiprocessing
savedir = const.txtdir
efficiency=np.loadtxt(savedir+'eff.txt')
a0=np.loadtxt(savedir+'a0.txt')
print('a0.argmax():',a0.argmax())
print('eff.argmax():',efficiency.argmax())
start = const.start
stop  = const.stop
step  = const.step
b = const.x_max/3e8/const.dt_snapshot/2
b = int(b)
print('b',b)
dirsdf  = const.sdfdir
dirsize =  const.filenumber
def rm(n):
	save = 0
	index_max = efficiency.argmax()+1
	a0_max    = a0.argmax()+1
	#print('a0',a0.max())
	#print('eff',efficiency.max())
	if n == b: 
		save = 1
	if n%100 == 0:
		save = 1
	if abs(n - index_max) % 50 == 0 and abs(n - index_max) < 2000:
                save = 1
	if abs(n - a0_max) % 50 == 0 and abs(n - a0_max) < 2000:
		save = 1
	if save == 0:
		if os.path.exists(dirsdf+str(n).zfill(dirsize)+".sdf")==True:
			print('begin_remove')
			os.remove(dirsdf+str(n).zfill(dirsize)+".sdf")
			print('remove:',n)
	if save == 1:
		print('save_sdf:',n)

	return
pool = multiprocessing.Pool(processes=4)
results=pool.map(rm,range(start,stop+step,step))
#print(len(results))


