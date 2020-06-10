
for n in range(start,stop+step,step):
        #### header data ####
        data = sdf.read(dirsdf+str(n).zfill(dirsize)+".sdf",dict=True)
        header=data['Header']
        time=header['time']
        Ex=data["Electric Field/Ex"].data
        Ex_y0=Ex[...,int(const.Ny/2)]
        Ey=data["Electric Field/Ey"].data
        Ey_y0=Ey[...,int(const.Ny/2)]
        ne=data['Derived/Number_Density/electron1'].data
        ne_y0=ne[...,int(const.Ny/2)]

	a=[Ey_y0,time]
        All_Ey_y0.append(a)
np.savetxt("./txt/"+const.data_name+"allEy_y0.txt",All_Ey_y0)
