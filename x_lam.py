# -- coding: utf-8 --
import numpy as np
import matplotlib.pyplot as plt
import sdf
import matplotlib.pyplot as pl
from matplotlib.ticker import MultipleLocator, FuncFormatter
from scipy.interpolate import interp1d
import constant as const
plt.switch_backend('agg')

xf=np.loadtxt(const.txtdir+'xf.txt')
#name = "/acceleration/"
###
locate  =  23000         #micron
#constant
c       =  3e8
micron  =  1e-6
lamada  =  const.lamada     #10.6 * micron
gridnumber = 2400
stop    =  const.stop     #21667
dt_snapshot=const.dt_snapshot #3e-15
dt      =  dt_snapshot*1e15      #fs
x_max   =  const.x_max      #60 * lamada
x_min   =  0 * lamada
x_end   =  x_max - x_min
window_start_time =  (x_max - x_min) / c
delta_x =  x_end/gridnumber
t_end   =  stop * dt_snapshot
x_interval=1
t_total=1e15*x_end/c         #fs
t_size=t_total/(dt_snapshot*1e15)+1+1           #t_grid_number
######t_size=int(1e15*gridnumber*delta_x/c)+1
x       = int(locate/(delta_x*x_interval*1e6))
#######
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
a=float("inf")
freqs[0]=a
light=3e8*np.ones(freqs.shape)
lam=(light/(freqs*1e12))*1e6
#lamda[1]=0
lam=lam[1:len(lam)-1]
#lamda[2]=0
#lamda[3]=0
sequence=len(lam)-1
#for i in range(len(lam)):
#   if lam[i]>10:
#      lam[i]=0
#for i in range(len(freqs)):
 #  if freqs[i]>500:
  #    freqs[i]=0
#####


#####set x ,y         


###interpolotion

####transition Xf
Xf=xf[x]
x_f=Xf[1:len(Xf)-1]
####
xx = np.linspace(lam.min(), lam.max(), 7000)
f = interp1d(lam, x_f,kind="cubic")
ynew=f(xx)
#plot

fig,ax=plt.subplots()
line=ax.plot(xx,ynew,"g")
#line2=ax.scatter(lam,x_f)

          
ax.set_xlabel('um')
ax.set_ylabel('')
plt.xlim((0,200))

#print and save
fig.savefig(const.figdir +  const.name +"_" +str(locate)+"lamda.png",dpi=200)
