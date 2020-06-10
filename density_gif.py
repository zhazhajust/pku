import numpy as np
import matplotlib.pyplot as plt
import sdf
import imageio
import os
import scipy.signal as signal
import constant as const
plt.switch_backend('agg')
###
start=1
stop =17000
step =1000
####
sdfdir="../Data/test_niezan/"
filenumber=5
savedir ="./gif/"
density_name="test_niezan.gif"
if (os.path.isdir(savedir) == False):
    os.mkdir(savedir)
image_list=[]
for i in range(start,stop,step):
  x=i
  savefigdir="./gif/png/"+str(x)
  data=sdf.read(sdfdir+str(x).zfill(filenumber)+".sdf",dict=True)
  Bz=data['Electric Field/Ey']
  time=data['Header']['time']
  bz=Bz.data
  density=data['Derived/Number_Density/electron1'].data
  bz=bz.T
  density=density.T
  bz_y0=bz[500]
  k0=2*3.14/0.8e-6
  f,t,zxx=signal.stft(bz_y0,fs=2*3.14/const.delta_x/k0,nperseg=const.nperseg)
  fig,axs=plt.subplots(3,1)
  im=axs[0].pcolormesh(bz,cmap=plt.get_cmap('bwr'))
  im2=axs[1].pcolormesh(density,cmap=plt.get_cmap('gray'))
  im3=axs[2].pcolormesh(t,f,np.abs(zxx),cmap=plt.get_cmap('BuPu'))
  axs[2].set_ylim((0,2))
  axs[0].set_title(time,fontsize=12,color='r')
  fig.savefig(savefigdir+"bz.png",dpi=200)
  
  image_list.append(str(savefigdir+"bz.png"))
  plt.cla()
  plt.clf()
  plt.close('all')
###
def create_gif(image_list, gif_name, duration = 1.0):
    '''
    ''' 
    frames = []
    for image_name in image_list:
        frames.append(imageio.imread(image_name))

    imageio.mimsave(gif_name, frames, 'GIF', duration=duration)
    return

def main():
    #image_list = ['1.jpg', '2.jpg', '3.jpg']
    gif_name = savedir + density_name
    duration = 1
    create_gif(image_list, gif_name, duration)

#if __name__ == '__main__':
main()
