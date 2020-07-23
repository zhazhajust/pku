def a0(x):
	dfdir=const.sdfdir +str(x).zfill(const.filenumber)+".sdf"
        data=sdf.read(sdfdir,dict=True)
        Ey=data['Electric Field/Ey']
        time=data['Header']['time']
        ey=Ey.data
	E0=ey.max()
	:
