import numpy as np
import constant as const
savedir = const.txtdir
eff=np.loadtxt(savedir+'eff.txt')
a0=np.loadtxt(savedir+'a0.txt')
a0_3 = np.loadtxt(savedir+'3.5Thz_a0.txt')
time1 =  np.loadtxt(const.txtdir + 'a0.txt').argmax()+1   #sdf_locate
locate1 = (time1*const.dt_snapshot - const.window_start_time)*3e8*1e6
###
time2 = np.loadtxt(const.txtdir + 'eff.txt').argmax()+1
locate2 = (time2*const.dt_snapshot - const.window_start_time)*3e8*1e6
print('a0.argmax():',time1,locate1,a0.max())
print('eff.argmax():',time2,locate2,eff.max())
print('a0_3.argmax():',a0_3.argmax(),(a0_3.argmax()*const.dt_snapshot - const.window_start_time)*3e8*1e6)
