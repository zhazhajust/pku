# -- coding: utf-8 --
import numpy as np
import matplotlib.pyplot as plt
import sdf
import matplotlib.pyplot as pl
from matplotlib.ticker import MultipleLocator, FuncFormatter
from scipy.interpolate import interp1d
import constant as const
plt.switch_backend('agg')

#xf=np.loadtxt(const.txtdir+'xf.txt')
#name = "/acceleration/"
###
locate  =  3400        #micron
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
eff = np.loadtxt(const.txtdir + 'eff.txt')
#time = np.arange(3116,3349)
locate = np.loadtxt(const.txtdir + 'eff_locate.txt')

#locate = (time*const.dt_snapshot - const.window_start_time)*3e8*1e6

#locate = int(locate)
#THz=[]
fig,axs =plt.subplots()
line=axs.plot(locate,eff,label='eff')
#ax2=ax.twinx()
#line2=ax.scatter(lam,x_f)
#ax.legend(loc='best')
#ax2.legend(loc='best')
#ax3.legend(loc='best')         
#axs[0].set_xlabel('density')
axs.set_ylabel('efficiency')
axs.set_xlabel('um')
axs.set_ylim([0,0.2])
#ax.set_ylabel('')
#plt.xlim((0,200))

#print and save
print(str(const.figdir +  const.name) +"_eff.png")
fig.savefig(const.figdir +  const.name +"_eff.png",dpi=400)
