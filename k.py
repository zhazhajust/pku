import numpy as np
import matplotlib.pyplot as plt
import sdf
import os
import scipy.signal as signal
import constant as const
import imageio
plt.switch_backend('agg')
###
####
interval=200
image_list=[]
def draw(x):
	savefigdir=const.figdir+str(x)
	data=sdf.read(const.sdfdir+str(x).zfill(const.filenumber)+".sdf",dict=True)
	Bz=data['Electric Field/Ey']
	time=data['Header']['time']
	bz=Bz.data
	density=data['Derived/Number_Density/electron1'].data
	bz=bz.T
	density=density.T
	index=int(const.Ny/2)
	bz_y0=bz[index]
	k0=2*3.14/const.lamada
	f,t,zxx=signal.stft(bz_y0,fs=2*3.14/const.delta_x/k0,nperseg=const.nperseg)
	fig,axs=plt.subplots(3,1)
	im=axs[0].pcolormesh(bz,cmap=plt.get_cmap('bwr'))
	im2=axs[1].pcolormesh(density,cmap=plt.get_cmap('gray'))
	im3=axs[2].pcolormesh(t,f,np.abs(zxx),cmap=plt.get_cmap('BuPu'))
	axs[2].set_ylim((0,2))
	axs[0].set_title(time,fontsize=12,color='r')
	fig.savefig(savefigdir+"_k.png",dpi=200)
	image_list.append(savefigdir+"_k.png")
	plt.close('all')
         
for i in range(1,const.stop,interval):
        draw(i)
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
    gif_name = "gif/"+ "k.gif"
    duration = 1
    create_gif(image_list, gif_name, duration)

#if __name__ == '__main__':
main()

