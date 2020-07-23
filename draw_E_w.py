import numpy as np
import matplotlib.pyplot as plt
import sdf
import os
import scipy.signal as signal
import constant as const
import imageio
import function as func
import scipy.fftpack as fftpack
from matplotlib.ticker import MultipleLocator, FuncFormatter
from multiprocessing.dummy import Pool as ThreadPool
#import wigner
plt.switch_backend('agg')
###
####
interval=100
image_list=[]
png_savedir = "./gif/png/"

def x_formatter(x, pos):
        a=(const.delta_x*x + const.c*T)   *1e6 
        return  "%d"%int(a)
def x_formatter2(x, pos):
	#locate = x * const.nperseg
        a=x + const.c*T * 1e6
        return  "%d"%int(a)
def draw(x):
	#print "draw",x
	savefigdir=const.figdir+str(x)
	sdfdir=const.sdfdir +str(x).zfill(const.filenumber)+".sdf"
	data=sdf.read(sdfdir,dict=True)
	#data=sdf.read(const.sdfdir+str(x).zfill(const.filenumber)+".sdf",dict=True)
	Bz=data['Electric Field/Ey']
	global T
	time=data['Header']['time']
	if time-const.window_start_time<0:
		T=0	
	else:
		T=time-const.window_start_time	
	bz=Bz.data
	density=data['Derived/Number_Density/electron1'].data
	bz=bz.T
	density=density.T
	index=int(const.Ny/2)
	bz_y0=bz[index]
	k0=2*3.14/const.lamada
	k1=2*3.14/1e-6
	fs=2*3.14/const.delta_x/k0
	f,t,zxx=signal.stft(bz_y0,fs=2*3.14/const.delta_x,nperseg=const.nperseg)
	print t,t.shape
	ne=data['Derived/Number_Density/electron1'].data
	ne_y0=ne[...,int(const.Ny/2)]
	fig,axs=plt.subplots(2,2)
	#plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

	ref=func.ref_a(x)
	ref=np.array(ref)
	im=axs[0][0].pcolormesh(bz,cmap=plt.get_cmap('bwr'))
	#fig.colorbar(im,ax=axs[0][0])
	#print ref
	line2=axs[1][0].plot(ref)     
        #(density,cmap=plt.get_cmap('gray'))
	im3=axs[0][1].pcolormesh(t*2*3.14/1e-6,f/k0,np.abs(zxx),cmap=plt.get_cmap('BuPu'),shading='gouraud')
	#Xf_list=wigner.wigner(x)

	#im3=axs[1][0].pcolormesh(Xf_list,cmap=plt.get_cmap('BuPu'),rasterized=True)
	#im=axs[1][0].imshow(np.abs(zxx),extent=[],cmap=plt.get_cmap('BuPu'))
	line3=axs[1][1].plot(ne_y0,'g')
	ax2=axs[1][1].twinx()
	#ax3=axs[0][1].twinx()
	hx = fftpack.hilbert(bz_y0)
	hy = np.sqrt(bz_y0**2+hx**2)
	line4=ax2.plot(hy,'r')
	#line5=axs[0][1].plot(ref)
	#ax3.set_ylim((0.9,1,1))
	
	axs[1][0].set_ylim((0.6,1.2))
	axs[1][0].set_ylabel('refractive_index')
	axs[0][1].set_ylim((0,2))
	axs[0][1].set_ylabel('k/k0')
	axs[1][1].set_ylabel('electron density')
	axs[1][1].set_xlabel("um")
	ax2.set_ylabel('electric field')
	axs[0][0].set_title("t="+str(time),fontsize=12,color='r')

	axs[0][0].xaxis.set_major_formatter( FuncFormatter( x_formatter ) )
	axs[0][1].xaxis.set_major_formatter( FuncFormatter( x_formatter2 ) )
	axs[1][0].xaxis.set_major_formatter( FuncFormatter( x_formatter ) )
	axs[1][1].xaxis.set_major_formatter( FuncFormatter( x_formatter ) )
	axs[1][1].xaxis.set_minor_locator( MultipleLocator(100) )
###
	fig.savefig(png_savedir+str(x)+"ref_k.png",dpi=200)
	image_list.append(png_savedir+str(x)+"ref_k.png")
	plt.close('all')


#wigner.aaa()         
for i in range(1,const.stop,interval):
        draw(i)



def create_gif(image_list, gif_name, duration = 0.5):
    '''
    '''
    frames = []
    #for image_name in image_list:
        #frames.append(imageio.imread(image_name))

    #imageio.mimsave(gif_name, frames, 'GIF', duration=duration)
 
    pool = ThreadPool()
    frames = pool.map(imageio.imread,image_list)
    imageio.mimsave(gif_name, frames, 'GIF', duration=duration)
    pool.close()
    pool.join()
    return

def main():
    #image_list = ['1.jpg', '2.jpg', '3.jpg']
    gif_name = "gif/"+const.name+".gif"
    duration = 0.5
    create_gif(image_list, gif_name, duration)

#if __name__ == '__main__':
const.checkdir()
main()
