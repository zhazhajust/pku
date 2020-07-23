# -- coding: utf-8 --
import numpy as np
import matplotlib.pyplot as plt
import sdf
import matplotlib.pyplot as pl
from matplotlib.ticker import MultipleLocator, FuncFormatter
import constant as const
plt.switch_backend('agg')

#name="contrast/density_half"
xf=np.loadtxt(const.txtdir+"xf.txt")
###
locate  =  22930        #micron

savedir = const.figdir  +  "freqs"+str(locate)+".png"
#constant
c       =  3e8
micron  =  1e-6
lamada  =  const.lamada   #10.6 * micron
gridnumber = 2400
stop    =  const.stop
dt_snapshot= const.dt_snapshot       #3e-15
dt      =  dt_snapshot*1e15      #fs
x_max   =  const.x_max        #100 * lamada
x_min   =  0 * lamada
x_end   =  x_max -x_min
window_start_time =  (x_max - x_min) / c
delta_x =  x_end/gridnumber
t_end   =  stop * dt_snapshot
x_interval=1
t_total=1e15*x_end/c         #fs
t_size=t_total/(dt_snapshot*1e15)+1+1           #t_grid_number
######t_size=int(1e15*gridnumber*delta_x/c)+1
x       = int(locate/(delta_x*x_interval*1e6))
print "x:"+str(x)
#######
if t_end-window_start_time<0:
      xgrid   =  int(gridnumber)
else:
      xgrid   =  int(gridnumber + c*(t_end-window_start_time)/delta_x)
#####fft freqs

N0 = t_size
T=t_size*dt             #fs  #dt_snapshot*1e15  #t[x][t_size-1]-t[x][0]
fs=N0*1e3/T
length=xf.shape[1]
freqs=np.linspace(0,fs/2,length)
################freqs=np.linspace(0,500,101)
#####time profile


#####


#####set x ,y         

####transition Xf
Xf=xf[x]
#plot
#np.savetxt(str(name)+str(locate)+"_freqs.txt",Xf)
fig,ax=plt.subplots()
line=ax.plot(freqs,Xf)
plt.xlim((0,10))
ax.xaxis.set_minor_locator( MultipleLocator(1) )
print "max:",str(max(Xf))         
ax.set_xlabel('Thz')
ax.set_ylabel('')


#print and save
fig.savefig(savedir,dpi=200)
