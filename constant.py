# -- coding: utf-8 --
import os
#constant
nperseg = 256
name    = 'a2_n3'
data_name = "a2_n3/"
filenumber = 4
sdfdir  =  "../Data/"+data_name
txtdir  =  "txt/"+data_name
figdir  =  "fig/"+data_name
def checkdir():
	if (os.path.isdir(txtdir) == False):
    		os.mkdir(txtdir)
	if (os.path.isdir(figdir) == False):
    		os.mkdir(figdir)
checkdir()
###
c       =  3e8
micron  =  1e-6
lamada  =  10.6 * micron 
gridnumber = 2400
Ny      =  2000
Nx      =  gridnumber
start   =  7000
stop    =  10000
step    =  1
dt_snapshot= 10e-15
dt      =  dt_snapshot*1e15      #fs
x_max   =  80 * lamada#* micron   #60 * lamada #micron
x_min   =  0  * lamada#* micron
x_end   =  x_max - x_min 
y_lenth =  100 * lamada
window_start_time =  (x_max - x_min) / c
delta_x =  x_end/gridnumber
t_end   =  stop * dt_snapshot
x_interval=1
t_total=1e15*x_end/c         #fs
t_size=t_total/(dt_snapshot*1e15)+1+1           #t_grid_number
######t_size=int(1e15*gridnumber*delta_x/c)+1

if t_end-window_start_time<0:
      xgrid   =  int(gridnumber)
else:
      xgrid   =  int(gridnumber + c*(t_end-window_start_time)/delta_x)
#####fft freqs
