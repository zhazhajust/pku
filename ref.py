import function as func 
import constant as const
import matplotlib.pyplot as plt 
from matplotlib.ticker import EngFormatter
from matplotlib.ticker import MultipleLocator, FuncFormatter

plt.switch_backend('agg')
savedir="./fig/"

t=2
def draw():
	a=func.ref_a(t)
	fig,ax=plt.subplots()
	###
	###
	line=ax.plot(a)
        ax.set_ylim((0.99,1))
	#fig.savefig('Xf.png',dpi=200)
	#set ticker

	def x_formatter(x, pos):
		a=const.delta_x*x*1e6
		return  "%d"%int(a)
	def freqs_formatter(x, pos):

		return  "%d"%int(x)
	x_major_locator=int(const.xgrid/const.x_interval/5)

	x_minor_locator=int(const.xgrid/const.x_interval/10)

	#y_tick_pos  = np.linspace(0,40,1)
	#ax.set_yticks(y_tick_pos)
	######
	#ax.set_yscale("symlog",basey=2)
	#ax.set_ylim((0,50))
	ax.xaxis.set_major_locator( MultipleLocator(x_major_locator) )
	ax.xaxis.set_major_formatter( FuncFormatter( x_formatter ) )
	ax.xaxis.set_minor_locator( MultipleLocator(x_minor_locator) )
	fig.savefig(savedir+const.data_name+str(t)+"ref.png",dpi=200)


draw()
