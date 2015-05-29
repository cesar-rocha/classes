import numpy as np
import matplotlib.pyplot as plt
import scipy.io as io

path = 'mae214A_SPR15_hw4/'
reystress = np.loadtxt(path+'chan395.reystress')
means = np.loadtxt(path+'chan395.means')

# reystress
# y,  y+ , R_uu, R_vv, R_ww, R_uv, R_uw , R_vw

# means
# y, y+, Umean, dUmean/dy, Wmean, dWmean/dy, Pmean
