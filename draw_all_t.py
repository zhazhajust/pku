# -- coding: utf-8 --
import numpy as np
import matplotlib.pyplot as plt
import sdf
import matplotlib.pyplot as pl
from matplotlib.ticker import MultipleLocator, FuncFormatter
plt.switch_backend('agg')
xt=np.loadtxt('./txt1/xt.txt')
npsave="./fig1/x_time.png"

#constant
c       =  3e8
micron  =  1e-6
gridnumber = 2400
stop    =  2000
dt_snapshot= 0.3e-15
dt      =  dt_snapshot*1e15      #fs
x_end   =  60 * micron
x_max   =  60 * micron
x_min   =  0 * micron
window_start_time =  (x_max - x_min) / c
delta_x =  x_end/gridnumber
t_end   =  stop * dt_snapshot
x_interval=1
t_total=1e15*x_end/c         #fs
t_size=t_total/(dt_snapshot*1e15)+1           #t_grid_number
######t_size=int(1e15*gridnumber*delta_x/c)+1

if t_end-window_start_time<0:
      xgrid   =  int(gridnumber)
else:
      xgrid   =  int(gridnumber + c*(t_end-window_start_time)/delta_x)
#####fft freqs

N0 = t_size
T=t_size*dt             #fs  #dt_snapshot*1e15  #t[x][t_size-1]-t[x][0]
fs=N0*1e3/T
freqs=np.linspace(0,fs/2,int(N0/2)+1)
################freqs=np.linspace(0,500,101)
#####time profile
t=np.arange(0,t_size+1,1)


#####


#####set x ,t        
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
fig.savefig(npsave,dpi=200)
