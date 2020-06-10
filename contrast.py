# -- coding: utf-8 --
import numpy as np
import matplotlib.pyplot as plt
import sdf
import matplotlib.pyplot as pl
from matplotlib.ticker import MultipleLocator, FuncFormatter
plt.switch_backend('agg')

xf=np.loadtxt('txt/density2.0/xf.txt')
xf2=np.loadtxt('txt/density_100_half/xf.txt')
###
locate  =  12000        #micron
savedir = "fig/contrast/freqs"+str(locate)+".png"
lenth   =  [100,60]         #box_lenth  /micron
#constant
c       =  3e8
micron  =  1e-6
lamada  =  10.6 * micron
gridnumber = 2400
stop    =  22967
dt_snapshot= 3e-15
dt      =  dt_snapshot*1e15 
for i in range(0,len(lenth))     #fs
   x_max[i]  = length[i] * lamada
   x_min   =  0 * lamada
   x_end[i]  =  x_max[i] -x_min
   window_start_time[i] =  (x_max[i] - x_min) / c
   delta_x[i] =  x_end[i]/gridnumber
   t_end   =  stop * dt_snapshot
   x_interval=10
   t_total[i]=1e15*x_end[i]/c         #fs
   t_size[i]=t_total[i]/(dt_snapshot*1e15)+1+1           #t_grid_number
######t_size=int(1e15*gridnumber*delta_x/c)+1
   x[i]= int(locate/(delta_x[i]*x_interval*1e6))
#######
if t_end-window_start_time[i]<0:
      xgrid   =  int(gridnumber)
else:
      xgrid   =  int(gridnumber + c*(t_end-window_start_time[i])/delta_x[i])
#####fft freqs

N0[i] = t_size[i]
T[i]=t_size[i]*dt             #fs  #dt_snapshot*1e15  #t[x][t_size-1]-t[x][0]
fs[i]=N0[i]*1e3/T[i]
freqs_length=xf[i].shape[1]
freqs=np.linspace(0,fs/2,freqs_length)
################freqs=np.linspace(0,500,101)
#####time profile


#####


#####set x ,y         

####transition Xf
Xf[i]=xf[x[i]]

freqs2=freqs
#plot
fig,ax=plt.subplots()
line=ax.plot(freqs,Xf)
line2=ax.plot(freqs2,Xf2,'r')
plt.xlim((0,50))
          
ax.set_xlabel('Thz')
ax.set_ylabel('')
ax.xaxis.set_minor_locator( MultipleLocator(2) )

#print and save
fig.savefig(savedir,dpi=200)
