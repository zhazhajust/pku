# -- coding: utf-8 --
import numpy as np
import matplotlib.pyplot as plt
import sdf
import constant as const
import matplotlib.pyplot as pl
from matplotlib.ticker import MultipleLocator, FuncFormatter
plt.switch_backend('agg')
###   
time =  np.loadtxt(const.txtdir + 'a0.txt').argmax()   #sdf_locate
locate = (time*const.dt_snapshot - const.window_start_time)*3e8*1e6
### 
time = np.loadtxt(const.txtdir + 'eff_Thz.txt').argmax()
locate2 = (time*const.dt_snapshot - const.window_start_time)*3e8*1e6
### 
### um  ########
x_locate=23000

########
###
locate='contrast' + str(x_locate)
xf=np.loadtxt(const.txtdir + 'xf.txt')

#xf=np.loadtxt('txt/density1e-2/xf.txt')
#xf2=np.loadtxt('txt/ac_a0_1/xf.txt')
#xf3=np.loadtxt('txt/density1.2/xf.txt')

savedir = const.figdir +locate +".png"
def xf_index(x_locate,x_lenth=80):
###
	locate  =  x_locate        #micron
	#savedir = "fig/density_100_half/freqs"+str(locate)+".png"
	#constant
	c       =  3e8
	micron  =  1e-6
	lamada  =  const.lamada #10.6 * micron
	gridnumber = 2400
	stop    =  const.stop #17667
	dt_snapshot= const.dt_snapshot   #3e-15
	dt      =  dt_snapshot*1e15      #fs

	x_max   =  x_lenth * lamada
	x_min   =  0 * lamada
	x_end   =  x_max -x_min
	window_start_time =  (x_max - x_min) / c
	delta_x =  x_end/gridnumber
	t_end   =  stop * dt_snapshot
	x_interval=const.x_interval
	t_total=1e15*x_end/c         #fs
	t_size=t_total/(dt_snapshot*1e15)+1+1           #t_grid_number
	######t_size=int(1e15*gridnumber*delta_x/c)+1
        x       = int(locate/(delta_x*x_interval*1e6))
        if t_end-window_start_time<0:
          xgrid   =  int(gridnumber)
        else:
          xgrid   =  int(gridnumber + c*(t_end-window_start_time)/delta_x)
	#####fft freqs

	N0 = t_size
	T=t_size*dt             #fs  #dt_snapshot*1e15  #t[x][t_size-1]-t[x][0]
	fs=N0*1e3/T
	#length=xf.shape[1]
	length=N0/2+1
        freqs=np.linspace(0,fs/2,length)
	return [x,freqs]
################freqs=np.linspace(0,500,101)
#####time profile



#####set x ,y         

####transition Xf
x=xf_index(x_locate)[0]
freqs=xf_index(x_locate)[1]
Xf=xf[x]
#plot
fig,ax=plt.subplots()
line=ax.plot(freqs,Xf,label='x1')
#plt.xlim((0,50))
#print "max:",str(max(Xf))         
ax.set_xlabel('Thz')
ax.set_ylabel('')
###
x=xf_index(x_locate+400)[0]
freqs2=xf_index(x_locate+400)[1] 
Xf2=xf[x]
line2=ax.plot(freqs2,Xf2,'g',label='x2')
###
x=xf_index(x_locate+800)[0]
freqs3=xf_index(x_locate+800)[1]
Xf3=xf[x]
line3=ax.plot(freqs3,Xf3,'g',label='x3')
ax.legend(loc='best')
ax.xaxis.set_minor_locator( MultipleLocator(1) )
plt.xlim((0,50))
#print and save
fig.savefig(savedir,dpi=200)
