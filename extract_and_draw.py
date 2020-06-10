# -- coding: utf-8 --
import math
import sdf
import matplotlib
import math
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from scipy.interpolate import spline
plt.switch_backend('agg')

xt_savedir="./txt/acceleration/"
lamadadir ="./fig/acceleration/"
xt_savename="xt.txt"
xf_savedir =  "./txt/acceleraton/xf.txt" 
lamada_savename="lamada.png"
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
  wavelength=     1.0e-6
  frequency =     v0*2*pi/wavelength
  micron    =     1.0e-6
  c         =     3e8
  exunit    =     m0*v0*frequency/q0
  bxunit    =     m0*frequency/q0
  denunit    =     frequency**2*epsilon0*m0/q0**2
  print 'electric field unit: '+str(exunit)
  print 'magnetic field unit: '+str(bxunit)
  print 'density unit nc: '+str(denunit)
  
  font = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 28,  
          }  
  
  
  
  ######### Parameter you should set ###########
  dirsdf  =  '../Data/acceleration/'
  dirsize =  5

#####
  micron  =  1e-6
  lamada  =  10.6 * micron
  c       =  3e8
  x_max   =  60 * lamada
  x_min   =  0 * micron
  window_start_time =  (x_max - x_min) / c 
  dt_snapshot= 3e-15  #fs
  dt      =  dt_snapshot*1e15  #fs
  start_move_number =  int(window_start_time / dt_snapshot)
  y       =  1250
  x_end   =  x_max - x_min
  gridnumber = 2400
  delta_x =  x_end/gridnumber
  start   =  1  # start time
  stop    =  22967 #21667  # end time
  step    =  1  # the interval or step
  t_end   =  stop * dt_snapshot
  t_n     =  int(t_end/1e-15)
  if t_end-window_start_time<0:
      xgrid   =  int(gridnumber)
  else:
      xgrid   =  int(gridnumber + c*(t_end-window_start_time)/delta_x)
   
  if (os.path.isdir(savedir) == False):
    os.mkdir(savedir)
  if (os.path.isdir(lamadadir) == False):
    os.mkdir(lamadadir)
  ######### Script code drawing figure ################
 
####################
  savedir
  x_interval=10
  t_total=1e15*x_end/c         #fs
  t_size=int(t_total/dt)+1+1   

######allay define
  xt=np.zeros((int(xgrid/x_interval)+1,t_size))

  for n in range(start,stop+step,step):
        #### header data ####
        data = sdf.read(dirsdf+str(n).zfill(dirsize)+".sdf",dict=True)
        header=data['Header']
        time=header['time']
        if  n  <  start_move_number:
                     
           for x in range(1,int(gridnumber/x_interval)+1):
              a=int(x*x_interval)
              d_n=int((1e15*delta_x*a/c)/dt)
              if n-d_n > 0 and n-d_n < t_size :
                    #[fs]
                   xt[x][n-d_n]=data['Magnetic Field/Bz'].data[a-1][y]/bxunit            
        else:
           for x in range(1,int(xgrid/x_interval)+1):
                
             #if x > 1200 :
               a=int(x*x_interval)
	       if a-c*(time-window_start_time)/delta_x >= 0 and a-c*(time-window_start_time)/delta_x < gridnumber-1:
		    #[fs]
                   d_n=int((1e15*delta_x*a/c)/dt)
                   xt[x][n-d_n]=data['Magnetic Field/Bz'].data[int(round(a-c*(time-window_start_time)/delta_x))][y]/bxunit
                   #else:bz.append(0)
                   #print 'Reading finished%d' %len(t)
  
  np.savetxt(xt_savedir+xt_savename, xt)
#############

xt=np.loadtxt(xt_savedir+xt_savename)

N0 = t_size
T=t_size*dt             #fs  #dt_snapshot*1e15  #t[x][t_size-1]-t[x][0]
fs=N0*1e3/T
freqs=np.linspace(0,fs/2,int(N0/2)+1)
################freqs=np.linspace(0,500,101)
#####time profile
t=np.arange(0,t_size+dt,dt)

rfft=np.fft.rfft(xt)
rfft=np.abs(rfft)
print("writed")
np.savetxt(xf_savedir, rfft)
print("saved")
xf=np.loadtxt('txt/xf.txt')
#########


x=np.arange(int(xgrid/x_interval)+1)
X,time=np.meshgrid(x,t)


####transition Xf
Xt=xt.T
#plot
fig,ax=plt.subplots()
im=ax.pcolormesh(Xt,cmap=plt.get_cmap('rainbow'))
fig.colorbar(im,ax=ax)
#fig.savefig('Xf.png',dpi=200)
#set ticker
def x_formatter(x, pos):
        delta_x=60*1e-6/2400
        a=delta_x*x*x_interval*1e6
        return  "%d"%int(a)
def t_formatter(x, pos):
          
        return  "%d"%int(x*dt)
x_major_locator=int(xgrid/x_interval/5)
x_minor_locator=int(xgrid/x_interval/10)
ax.xaxis.set_major_locator( MultipleLocator(x_major_locator) )
ax.xaxis.set_major_formatter( FuncFormatter( x_formatter ) )
ax.xaxis.set_minor_locator( MultipleLocator(x_minor_locator) )
ax.yaxis.set_major_formatter( FuncFormatter( t_formatter ) )
ax.set_xlabel('um')
ax.set_ylabel('fs')



#print and save
plt.show()
fig.savefig("fig/x_time.png",dpi=200)


########
#####set x ,y         
x=np.arange(int(xgrid/x_interval)+1)
a=float("inf")
freqs[0]=a
light=3e8*np.ones(freqs.shape)
lam=(light/(freqs*1e12))*1e6
X,lamda=np.meshgrid(x,lam)
lamda[1]=0
lamda[0]=0
lamda[2]=0
lamda[3]=0

####transition Xf
Xf=xf.T
#plot
fig,ax=plt.subplots()
im=ax.pcolormesh(X,lamda,Xf,cmap=plt.get_cmap('rainbow'),shading='gouraud')
#im=ax.imshow(X,lamda,Xf,cmap=plt.get_cmap('rainbow'),interpolation='sinc')
fig.colorbar(im,ax=ax)
#fig.savefig('Xf.png',dpi=200)
#set ticker
def x_formatter(x, pos):
          delta_x=60*1e-6/2400
          a=delta_x*x_interval*x*1e6
          return  "%d"%a

x_major_locator=int(xgrid/x_interval/5)
x_minor_locator=int(xgrid/x_interval/10)
ax.xaxis.set_major_locator( MultipleLocator(x_major_locator) )
ax.xaxis.set_major_formatter( FuncFormatter( x_formatter ) )
ax.xaxis.set_minor_locator( MultipleLocator(x_minor_locator) )
ax.set_xlabel('um')
ax.set_ylabel('um')
#print and save
plt.show()
fig.savefig(lamadadir+ lamada_savename,dpi=200)
