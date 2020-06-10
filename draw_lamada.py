# -- coding: utf-8 --
import numpy as np
import matplotlib.pyplot as plt
import sdf
import matplotlib.pyplot as pl
from matplotlib.ticker import MultipleLocator, FuncFormatter
plt.switch_backend('agg')

xf=np.loadtxt('txt/a0_1_n^2/xf.txt')
savedir = "fig/a0_1_n^2/xf.png"
#constant
c       =  3e8
micron  =  1e-6
lamada  =  10.6 * micron
gridnumber = 2400
stop    =  35334
dt_snapshot= 3e-15
dt      =  dt_snapshot*1e15      #fs

x_max   =  80 * lamada
x_min   =  0 * lamada
x_end   =  x_max-x_min
window_start_time =  (x_max - x_min) / c
delta_x =  x_end/gridnumber
t_end   =  stop * dt_snapshot
x_interval=10
t_total=1e15*x_end/c         #fs
t_size=t_total/(dt_snapshot*1e15)+1+1           #t_grid_number
######t_size=int(1e15*gridnumber*delta_x/c)+1

if t_end-window_start_time<0:
      xgrid   =  int(gridnumber)
else:
      xgrid   =  int(gridnumber + c*(t_end-window_start_time)/delta_x)
#####fft freqs

N0 = t_size
T=t_size*dt             #fs  #dt_snapshot*1e15  #t[x][t_size-1]-t[x][0]
fs=N0*1e3/T
lenth = xf.shape[1]
freqs=np.linspace(0,fs/2,N0/2+1)    #int(N0/2)+1)
################freqs=np.linspace(0,500,101)
#####time profile
t=np.arange(0,t_size+dt,dt)


#####


#####set x ,y         
x=np.arange(int(xgrid/x_interval)+1)
a=float("inf")
freqs[0]=a
light=3e8*np.ones(freqs.shape)
lam=(light/(freqs*1e12))*1e6
X,lamda=np.meshgrid(x,lam)
#lamda[1]=0
lamda[0]=0
#lamda[2]=0
#lamda[3]=0
#lamda[4]=0
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
          a=delta_x*x_interval*x*1e6
          return  "%d"%a

x_major_locator=int(xgrid/x_interval/5)
x_minor_locator=int(xgrid/x_interval/10)
ax.set_yscale("symlog")
ax.set_ylim((5,50))
ax.xaxis.set_major_locator( MultipleLocator(x_major_locator) )
ax.xaxis.set_major_formatter( FuncFormatter( x_formatter ) )
ax.xaxis.set_minor_locator( MultipleLocator(x_minor_locator) )
ax.set_xlabel('um')
ax.set_ylabel('um')
#print and save
plt.show()
fig.savefig(savedir,dpi=200)
