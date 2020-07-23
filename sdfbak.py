import numpy as np
import matplotlib.pyplot as plt
import sdf
import os
import scipy.signal as signal
import constant as const
import function as func
import scipy.fftpack as fftpack
from matplotlib.ticker import MultipleLocator, FuncFormatter
plt.switch_backend('agg')
limit_min=0.1e12
limit_max=10e12
locate=12500e-6
locate2=16000e-6
middle = (locate/3e8 + const.window_start_time)/const.dt_snapshot
middle = int(middle)
d_n = 800e-6/3e8/const.dt_snapshot
d_n = int(d_n)
a1 = (locate/3e8 + const.window_start_time)/const.dt_snapshot
a2 = ((locate2-800e-6)/3e8 + const.window_start_time)/const.dt_snapshot
b = const.x_max/3e8/const.dt_snapshot/2
b = int(b)
print 'b',b
print 'sdf',a1,a2
print 2*d_n
