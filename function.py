import numpy
import sdf
import constant as const
import scipy.signal as signal
import numpy as np
pi=3.14
def E_x_y_zxx(a):        
	sdfdir=const.sdfdir +str(a).zfill(const.filenumber)+".sdf"
	data=sdf.read(sdfdir,dict=True)
	Ex=data["Electric Field/Ex"].data
	Ex_y0=Ex[...,int(const.Ny/2)]
	Ey=data["Electric Field/Ey"].data
	Ey_y0=Ey[...,int(const.Ny/2)]
	k,x,zxx=signal.stft(Ey_y0,fs=2*pi/const.delta_x,nperseg=const.nperseg)
	zxx=abs(zxx)
	index = np.unravel_index(zxx.argmax(),zxx.shape)
	k_x=np.ones(const.Nx)*k[index[0]]
	for i in range(0,const.Nx):
		x=i
		a=zxx[...,int(x/const.nperseg)]
		a=a.tolist()
		a[0]=0
		max_index=a.index(max(a))
		if max(a) > 0.2 * zxx[index[0]][index[1]]:
			k_x[i]=(k[max_index])
	ne=data['Derived/Number_Density/electron1'].data
	ne_y0=ne[...,int(const.Ny/2)]

	return [Ex_y0,Ey_y0,zxx,ne_y0,k_x]
def k(x,zxx):
	a=zxx[...,int(x/const.nperseg)]
	a=a.tolist()
	a[0]=0
	max_index=a.index(max(a))
	# index=(max_index*pi)/(len(a)*const.delta_x)
	index =  (pi/const.delta_x)/len(a) * max_index
	#print "max(k_x),max_index,len(k_x),index",max(a),max_index,len(a),index # ,index
	return index
def scalar_p(Ex_y0):
	scalar_p=0
	e=1.6e-19
	c=3e8
	m0=9.1e-31
	import constant as const
	a=np.zeros(const.Nx)
	x=const.Nx - 1
	for i in range(0,x+1):
		scalar_p=Ex_y0[x-i]*const.delta_x+scalar_p 
	#a.append(e*scalar_p/(m0*c**2))   
		a[x-i]=e*scalar_p/(m0*c**2)  
	#   print  'n_scalar',a
	return a
def wp_2(x,ne_y0):
	ne=ne_y0[x]
	#ne=1e25
	e=1.6e-19
	m0=9.1e-31
#	print "ne",ne
	w_p_2=ne*e**2/(m0*8.85e-12)
#	print "w_p",w_p_2
	return w_p_2
def ref_index(x,Ex_y0,Ey_y0,zxx,ne_y0,k_x):
	c=3e8
	wp_x_2=wp_2(x,ne_y0)
#	print "wp^2",wp_x_2
	#p_y0=scalar_p(Ex_y0)
	p_y0_x=p_y0[x]
#	print "scalar",p_y0_x
	#k_x=k(x,zxx)
	k_x=k_x[x]
	w0=k_x*c/1
#	print "w0",w0
	a=wp_x_2/(1+p_y0_x)
	b=(k_x*c)**2
	c=-(k_x*c)**2
#	print "a,b,c",a,b,c
	return np.roots([a,b,c])
def ref_a(a):   
	ref=[]
	Ex_y0,Ey_y0,zxx,ne_y0,k_x=E_x_y_zxx(a)
	x=const.Nx
	global p_y0
	p_y0=scalar_p(Ex_y0)
	for b in range(0,x):
		c=ref_index(b,Ex_y0,Ey_y0,zxx,ne_y0,k_x)
		#print("c",c)
		a= type(c) == numpy.ndarray
		#print("type",a)
		if type(c) == numpy.ndarray:
			c=c[-1]
		#print("c[-1]",c)
		ref.append(c)
	#print("ref:",type(ref[0]))
	return ref
#def help():
#	print "func.ref_a(a)"
#	print "a,x=1,1"
#	print "Ex_y0,Ey_y0,zxx,ne_y0,k_x=func.E_x_y_zxx(a)"
#	print "global p_y0"
#	print "func.p_y0=scalar_p(Ex_y0)"    
#	print "func.ref_index(x,Ex_y0,Ey_y0,zxx,ne_y0)"
