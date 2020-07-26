# -- coding: utf-8 --
import math
import sdf
import matplotlib
import math
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import multiprocessing 
import time
from numpy import ma
from matplotlib import colors, ticker, cm
import constant as const
import os
#from multiprocessing.dummy import Pool as ThreadPool 
savedir=const.txtdir   #"./txt/a0_1_2e-2/"
savename="xt.txt"
fftdir =const.figdir      #"./fig/a0_1_2e-2/"  
###
dirsdf  = const.sdfdir   # '../Data/a0_1_2e-2/'
dirsize =  const.filenumber    #4

bakdir="baktxt/"+const.data_name

b = const.x_max/3e8/const.dt_snapshot/2
b = int(b)
if (os.path.isdir(bakdir) == False):
        os.mkdir(bakdir)

def rm(n):
	save = 0
	#index_max = efficiency.argmax()+1
	#a0_max    = a0.argmax()+1
	if n==b:
		save = 1
	if n%50 == 0:
		save = 1
	if abs(n - index_max) % 10 == 0 and abs(n - index_max) < 1000:
		save = 1
	if abs(n - a0_max) % 10 == 0 and abs(n - a0_max) < 1000:
		save = 1	
	if save == 0:
		if os.path.exists(dirsdf+str(n).zfill(dirsize)+".sdf")==True:
			print('begin_remove')
			os.remove(dirsdf+str(n).zfill(dirsize)+".sdf")
			print('remove:',i)
	'''
	if save == 1:
		print('save_sdf:',i)
	'''	
	return
def extract_a0(x):
	limit_min=1e12
	limit_max=5.5e12
	freqs = 3.5e12
	c=3e8
	#w0=4e12
	me=9.10956e-31
	e=1.602176565e-19

	#p "draw",x

	#sdfdir=const.sdfdir +str(x).zfill(const.filenumber)+".sdf"
	#data=sdf.read(sdfdir,dict=True)

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
	#print('E0:'+str(E0_w0))
	w0=freqs*2*3.1415926
	a0_w0=e*E0_w0/(me*c*w0)
	return a0_w0
def energe(x):
        #p "draw",x
        savefigdir=const.figdir+str(x)+'k_bz.png'

        #sdfdir=const.sdfdir +str(x).zfill(const.filenumber)+".sdf"
        #data=sdf.read(sdfdir,dict=True)

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
def baktxt(n):
        #print("save file:"+str(n))

        #sdfdir=const.sdfdir +str(n).zfill(const.filenumber)+".sdf"
        #data = sdf.read(sdfdir,dict=True)

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

def extract(n):
        #### header data ####
	print('n:'+str(n))
	#xt = np.ndarray(ss1, dtype=ss2, buffer=shm.buf)
	data = sdf.read(dirsdf+str(n).zfill(dirsize)+".sdf",dict=True)
	header=data['Header']
	time=header['time']
	E_y0=data['Electric Field/Ey'].data[:,int(y)]
	if  n  <  start_move_number:
		for x in range(1,int(gridnumber/x_interval)+1):
			a=int(x*x_interval)
			d_n=int((1e15*delta_x*a/c)/dt)
			if n-d_n > 0 and n-d_n < t_size :
#[fs]
				xt[x][n-d_n]=E_y0[a-1] #/bxunit            
	else:
		for x in range(1,int(xgrid/x_interval)+1):

			a=int(x*x_interval)
			if a-c*(time-window_start_time)/delta_x >= 0 and a-c*(time-window_start_time)/delta_x < gridnumber-1:
	#[fs]
				d_n=int((1e15*delta_x*a/c)/dt)
				xt[x][n-d_n]=E_y0[int(round(a-c*(time-window_start_time)/delta_x))]  #/bxunit
	return "OK"+str(n)
if __name__ == "__main__":
	######## Constant defined here ########
	pi        =     3.1415926535897932384626
	q0        =     1.602176565e-19 # C
	m0        =     9.10938291e-31  # kg
	v0        =     2.99792458e8    # m/s^2
	kb        =     1.3806488e-23   # J/K
	mu0       =     4.0e-7*pi       # N/A^2
	epsilon0  =     8.8541878176203899e-12 # F/m
	h_planck  =     6.62606957e-34  # J s
	####lamada


	wavelength=     const.lamada     #10.6e-6

	####

	frequency =     v0*2*pi/wavelength
	micron    =     1e-6
	c         =     3e8
	exunit    =     m0*v0*frequency/q0
	bxunit    =     m0*frequency/q0
	denunit    =     frequency**2*epsilon0*m0/q0**2
	print('electric field unit: '+str(exunit))
	print('magnetic field unit: '+str(bxunit))
	print('density unit nc: '+str(denunit))
	font = {'family' : 'monospace',  
		'color'  : 'black',  
		'weight' : 'normal',  
		'size'   : 28,  
	}  
	if (os.path.isdir(savedir) == False):
		os.mkdir(savedir)
		
	if (os.path.isdir(fftdir) == False):
		os.mkdir(fftdir)    
	######### Script code drawing figure ################
	######constant
	###
	c       =  3e8
	micron  =  1e-6
	lamada  =  const.lamada #10.6 * micron
	gridnumber = const.Nx     #2400
	start   =  const.start
	stop    =  const.stop       #5889 #17000
	step    =  1
	dt_snapshot= const.dt_snapshot     #9e-15
	dt      =  dt_snapshot*1e15      #fs
	x_max   =  const.x_max      #80 * lamada   #60 * lamada    #micron
	x_min   =  0 * micron
	x_end   =  x_max - x_min
	y       =  const.Ny/2 
	window_start_time =  (x_max - x_min) / c
	#start_move_number = window_start_time * 1e15      #fs
	start_move_number =  int(window_start_time / dt_snapshot)
	delta_x =  x_end/gridnumber
	t_end   =  stop * dt_snapshot
	x_interval=const.x_interval          #10
	t_total=1e15*x_end/c         #fs
	t_size=t_total/(dt_snapshot*1e15)+1           #t_grid_number
	if t_end-window_start_time<0:
		xgrid   =  int(gridnumber)
	else: 
		xgrid   =  int(gridnumber + c*(t_end-window_start_time)/delta_x)

####################
	x_interval=const.x_interval      #10
	t_total=1e15*x_end/c         #fs
	t_size=int(t_total/dt)+1+1   

	SHAPE = ((int(xgrid/x_interval)+1,t_size))
	xt=np.zeros((int(xgrid/x_interval)+1,t_size))

	limit_min=0.1e12
	limit_max=10e12
	b = const.x_max/3e8/const.dt_snapshot/2
	b = int(b)
	print('b',b)
	e_start=float("inf")
	efficiency = np.zeros(int(stop))
	a0 = np.zeros(int(stop))
	for n in range(start,stop+step,step):
		#print('n',n)
		while 1:
			if os.path.exists(dirsdf+str(n+1).zfill(dirsize)+".sdf") == True:# or os.path.exists(dirsdf+str(stop).zfill(dirsize)+".sdf") == True:
				#time.sleep()
				break
		sdfdir=const.sdfdir +str(n).zfill(const.filenumber)+".sdf"
		data=sdf.read(sdfdir,dict=True)
		extract(n)
		baktxt(n)
		a0[n-1] = extract_a0(n)
		#print(a0)
		if n==b:
			e_start=energe(b)[0]		
		n_eff=energe(n)[1]/e_start
		efficiency[n-1]=n_eff
		if n % 10 == 0:
			index_max = efficiency.argmax()+1
			a0_max    = a0.argmax()+1
			print(range(1,n+1)[-1])
			for i in range(1,n+1):
				rm(i)			
	np.savetxt(savedir+savename, xt)
	np.savetxt(savedir+'eff.txt', efficiency)
	np.savetxt(savedir+'3.5Thz_a0.txt', a0)

	print('finished')
	#np.savetxt(savedir+savename, xt)
	#np.savetxt(savedir+'eff.txt', efficiency)
	#np.savetxt(savedir+'3.5Thz_a0.txt', a0)


	#time = range(0,stop+step,step)
	#locate = (time*const.dt_snapshot - const.window_start_time)*3e8*1e6
	#np.savetxt(const.txtdir + 'eff_locate.txt',locate)
