import numpy as np
import matplotlib.pyplot as plt
import sdf
import imageio
import os
import scipy.signal as signal
import constant
plt.switch_backend('agg')
###
start=1
stop =5000
step =200
####
sdfdir="../Data/density2e-2/"
filenumber=4
savedir ="./gif/"
density_name="density2e-2.gif"
E_name="e_density2e-2.gif"
if (os.path.isdir(savedir) == False):
    os.mkdir(savedir)
image_list=[]
for i in range(start,stop,step):
  x=i
  savefigdir="./gif/png/"+str(x)
  data=sdf.read(sdfdir+str(x).zfill(filenumber)+".sdf",dict=True)
  Bz=data['Magnetic Field/Bz']
  time=data['Header']['time']
  bz=Bz.data
  density=data['Derived/Number_Density/electron1'].data
  bz=bz.T
  bz_y0=bz[1000]
  f,x,zxx=signal.stft(bz_y0,fs=delta_x/c)
  density=density.T
  fig,axs=plt.subplots(3,1)
  im=axs[0].pcolormesh(bz,cmap=plt.get_cmap('bwr'))
  im2=axs[1].pcolormesh(density,cmap=plt.get_cmap('gray'))
  im3=axs[2].pcolormesh(x,f,zxx,cmap=plt.get_cmap('bwr'))
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
