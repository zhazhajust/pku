import constant as const
import sdf
import os
import numpy as np
start=const.start
stop=const.stop
step=const.step

savedir="baktxt/"+const.data_name

if (os.path.isdir(savedir) == False):
	os.mkdir(savedir)

for n in range(start,stop+step):
	print "save file",n
        sdfdir=const.sdfdir +str(n).zfill(const.filenumber)+".sdf"
        data = sdf.read(sdfdir,dict=True)
        header=data['Header']
        time=header['time']	
        Ex=data["Electric Field/Ex_averaged"].data
        Ex_y0=Ex[...,int(const.Ny/2)]
        Ey=data["Electric Field/Ey_averaged"].data
        Ey_y0=Ey[...,int(const.Ny/2)]
        ne=data['Derived/Number_Density/electron1'].data
        ne_y0=ne[...,int(const.Ny/2)]
	np.savetxt("baktxt/"+const.data_name+str(n)+"Ey_y0.txt",Ey_y0)
	np.savetxt("baktxt/"+const.data_name+str(n)+"Ex_y0.txt",Ex_y0)
	np.savetxt("baktxt/"+const.data_name+str(n)+"ne_y0.txt",ne_y0)
