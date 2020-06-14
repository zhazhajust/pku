# -- coding: utf-8 --
import sdf
import numpy as np
import constant as const
###
txtdir=const.txtdir
xt=np.loadtxt(txtdir+"xt.txt")
savedir =  txtdir + "xf.txt"
###
c       =  3e8
micron  =  1e-6
lamada  =  const.lamada #10.6 * micron
gridnumber = const.Nx     #2400
stop    =  const.stop       #5889 #17000
dt_snapshot= const.dt_snapshot     #9e-15
dt      =  dt_snapshot*1e15      #fs
x_max   =  const.x_max      #80 * lamada   #60 * lamada    #micron
x_min   =  0 * micron
x_end   =  x_max - x_min
window_start_time =  (x_max - x_min) / c
delta_x =  x_end/gridnumber
t_end   =  stop * dt_snapshot
x_interval=const.x_interval          #10
t_total=1e15*x_end/c         #fs
t_size=t_total/(dt_snapshot*1e15)+1           #t_grid_number
if t_end-window_start_time<0:
      xgrid   =  int(gridnumber)
else:
      xgrid   =  int(gridnumber + c*(t_end-window_start_time)/delta_x)
#####fft freqs

N0 = t_size
T=t_size*dt             #fs  #dt_snapshot*1e15  #t[x][t_size-1]-t[x][0]
fs=N0*1e3/T
rfft=np.fft.rfft(xt)
rfft=np.abs(rfft)
#f_size=rfft.shape
#xf=np.zeros((int(xgrid/x_interval)+1,f_size))
#for x in range(1,xf.shape[0]):
#    xf[x] = np.fft.rfft(xt[x])/N0
    #xf[0] = xf[0]/2
    #xf[N0/2] = xf[N0/2]/2
#    xf=np.abs(xf)
print("writed")
np.savetxt(savedir, rfft)
print("saved")

