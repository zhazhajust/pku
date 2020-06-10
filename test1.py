import imageio
import constant as const
interval=200
image_list=[]
png_savedir=",/fig/png/"
for i in range(int(const.stop/interval)):
   x=i
   #savefigdir=const.figdir+str(x)
   image_list.append(png_savedir+"ref_k.png")
   #image_list.append(savefigdir+"_k.png")
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

main()
