# -- coding: utf-8 --
import numpy as np
import matplotlib.pyplot as plt
import sdf
import matplotlib.pyplot as pl
from matplotlib.ticker import EngFormatter
from matplotlib.ticker import MultipleLocator, FuncFormatter
import constant as const
plt.switch_backend('agg')
load_dir=const.txtdir+"xf.txt"
xf=np.loadtxt(load_dir)
savedir=const.figdir +  "freqs.png"
#constant
c       =  3e8
micron  =  1e-6
lamada  =  const.lamada #10.6 * micron 
gridnumber = const.Nx      #2400
stop    =  const.stop      #5889
dt_snapshot= const.dt_snapshot      #9e-15
dt      =  dt_snapshot*1e15      #fs
x_max   =  const.x_max        #80 * lamada   #60 * lamada #micron
x_min   =  0 * micron
x_end   =  x_max - x_min 
window_start_time =  (x_max - x_min) / c
delta_x =  x_end/gridnumber
t_end   =  stop * dt_snapshot
x_interval=const.x_interval        #10
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
freqs=np.linspace(0,fs/2,int(N0/2)+1)
print(len(freqs))
######
for i in range(0,len(freqs)):
     if freqs[i] > 50:
         index = i
         break;

#freqs=freqs[0:index]

################freqs=np.linspace(0,500,101)
#####time profile
t=np.arange(0,t_size+dt,dt)


#####


#####set x ,y         
print(xgrid*delta_x*1e6/x_interval,delta_x*1e6)
#x=np.arange(int(xgrid*delta_x*1e6/x_interval)+1,delta_x*1e6)
x=np.linspace(0,xgrid*delta_x*1e6/x_interval,int(xgrid/x_interval+1))

#a=float("inf")
#freqs[0]=a
#light=3e8*np.ones(freqs.shape)
#lam=(light/(freqs*1e12))*1e6
X,Freqs=np.meshgrid(x,freqs)
#lamda[1]=0
#lamda[0]=0
#lamda[2]=0
#lamda[3]=0

####transition Xf
Xf=xf.T
#plot
fig,ax=plt.subplots()
###
#Xf=Xf[0:index,...]
###
im=ax.pcolormesh(X,Freqs,Xf,cmap=plt.get_cmap('rainbow'))
fig.colorbar(im,ax=ax)
#fig.savefig('Xf.png',dpi=200)
#set ticker

def x_formatter(x, pos):
        a=delta_x*x*x_interval*1e6
        return  "%d"%int(a)
def freqs_formatter(x, pos):
          
        return  "%d"%int(x)
x_major_locator=int(xgrid/x_interval/5)
x_minor_locator=int(xgrid/x_interval/50)

#y_tick_pos  = np.linspace(0,40,1)
#ax.set_yticks(y_tick_pos)
######
#ax.set_yscale("symlog",basey=2)
#ax.set_ylim((0,50))
#ax.xaxis.set_major_locator( MultipleLocator(x_major_locator) )
#ax.xaxis.set_major_formatter( FuncFormatter( x_formatter ) )
#ax.xaxis.set_minor_locator( MultipleLocator(x_minor_locator) )
ax.xaxis.set_minor_locator( MultipleLocator(200) )
ax.yaxis.set_minor_locator( MultipleLocator(1) )
#formatter0 = EngFormatter(unit='THz')
#ax.yaxis.set_major_formatter(formatter0)
ax.yaxis.set_major_formatter( FuncFormatter( freqs_formatter ) )
ax.set_xlabel('um')
ax.set_ylabel('Thz')
#print and save
#plt.show()
fig.savefig(savedir,dpi=200)
