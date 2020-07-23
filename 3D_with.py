import os
os.path.exists(test_file.txt)
def baktxt(n):
        print("save file:"+n)
        sdfdir=const.sdfdir +str(n).zfill(const.filenumber)+".sdf"
        data = sdf.read(sdfdir,dict=True)
        header=data['Header']
        time=header['time']
        Ex=data["Electric Field/Ex"].data
        Ex_y0=Ex[...,int(const.Ny/2)]
        Ey=data["Electric Field/Ey"].data
        Ey_y0=Ey[...,int(const.Ny/2)]
        ne=data['Derived/Number_Density/electron1'].data
        ne_y0=ne[...,int(const.Ny/2)]
        np.savetxt("baktxt/"+const.data_name+str(n)+"Ey_y0.txt",Ey_y0)
        np.savetxt("baktxt/"+const.data_name+str(n)+"Ex_y0.txt",Ex_y0)
        np.savetxt("baktxt/"+const.data_name+str(n)+"ne_y0.txt",ne_y0)

energe_start = total_energe(141)
n_sdf = 1
end_sdf = 

def draw(x):
        savefigdir=const.figdir+str(x)+'k_bz.png'
        sdfdir=const.sdfdir +str(x).zfill(const.filenumber)+".sdf"
        data=sdf.read(sdfdir,dict=True)
        Bz=data['Electric Field/Ey']
        time=data['Header']['time']
        bz=Bz.data
        bz=bz.T
        k_bz=np.fft.fft(bz)
        delta_k=3.14/const.delta_x/(const.Nx/2)
        k_bz2=k_bz*1
        k_n=[]
        for n in range(0,const.Nx):
                mi = 3e8/limit_min
                ma = 3e8/limit_max
                if 2 * 3.14 / ma  > n * delta_k and  n * delta_k > 2 * 3.14 / mi:
                        k_n.append(n)
        k_bz2[...,0:k_n[0]]=0    #k_bz.argmin()
        k_bz2[...,k_n[-1]:-k_n[-1]]=0  #k_bz.argmin()
        k_bz2[...,-k_n[0]:]=0    #k_bz.argmin()
	bz_filter=np.fft.ifft(k_bz2)
        E_x=np.sum(np.sum(np.square(bz)))
        E_Thz=np.sum(np.sum(np.square(bz_filter.real)))
        eff=E_Thz/E_x
        return [E_x,E_Thz]

def total_energe(n_sdf):
	t_e = draw(n_sdf)
	e1=t_e[0]
	return e1
def exist(n_sdf):
	if os.path.exists(str(n_sdf).zfill(5)+'.sdf') == True :
		if n_sdf = end_sdf:
			state = 100
		else:
			state = 1
	return state
def energe(n_sdf):
	e1=draw(n_sdf)
	e2=e1[1]
	return e2
def bak_sdf(n_sdf):
	baktxt(n_sdf)

while 1:
	if exist(n_sdf) = 1:
		if txt_exist(n_sdf) = 0:
			bak_sdf(n_sdf)
			efficiency = energe(n_sdf,f=[0.1,10])/energe_start
			eff.append(efficiency)
			np.savetxt('eff'+str(n_sdf)+'.txt',eff)
			if exist('eff'+str(n_sdf-1)+'.txt'):
				rm('eff'+str(n_sdf-1)+'.txt')			
			if eff[n-1] != eff.max() and  n_sdf%100 != 0 :
				rm(n_sdf)
		else:
			n_sdf = n_sdf + 1
	if exist(n_sdf) = 100:
		break
	else:
		sleep(60)
		continue

			
				
