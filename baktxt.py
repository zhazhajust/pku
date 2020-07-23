import constant as const
import sdf
import os
import numpy as np
import multiprocessing
start=const.start
stop=const.stop
step=const.step

savedir="baktxt/"+const.data_name

if (os.path.isdir(savedir) == False):
	os.mkdir(savedir)
#for i in range(start,stop+step,step):
#       results.append(pool.apply_async(extract, (i, ))) 

'''
for n in range(start,stop+step):
	baktxt(n)
'''
def baktxt(n):
	print("save file:"+str(n))
	sdfdir=const.sdfdir +str(n).zfill(const.filenumber)+".sdf"
	data = sdf.read(sdfdir,dict=True)
	header=data['Header']
	time=header['time']	
	Ex=data["Electric Field/Ex"].data
	Ex_y0=Ex[...,int(const.Ny/2)]
	Ey=data["Electric Field/Ey"].data
	Ey_y0=Ey[...,int(const.Ny/2)]
	ne=data['Derived/Number_Density/electron1'].data
	ne_y0=ne[...,int(const.Ny/2)]
	np.savetxt("baktxt/"+const.data_name+str(n)+"Ey_y0.txt",Ey_y0)
	np.savetxt("baktxt/"+const.data_name+str(n)+"Ex_y0.txt",Ex_y0)
	np.savetxt("baktxt/"+const.data_name+str(n)+"ne_y0.txt",ne_y0)

pool = multiprocessing.Pool(processes=4)
#for i in range(start,stop+step,step):
#       results.append(pool.apply_async(extract, (i, ))) 
results = pool.map(baktxt,range(int(start),int(stop+step)))

pool.close()
pool.join()
