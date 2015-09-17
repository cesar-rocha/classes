
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from mpltools import style


plt.close('all')
plt.rcParams.update({'font.size': 30, 'legend.handlelength'  : 1.5
        , 'legend.markerscale': 14., 'legend.linewidth': 3.})

#mu2 = np.linspace(0.,2.,)
mu = 15.
z = np.linspace(0.,1.,30)

phi1 = sp.exp(-mu*z)
phi2 = sp.exp(mu*(z-1))
color1 = '#ffd700'
color3 = '#9acd32'
color2 = '#ff0000'

lw1=4
aph=.7

#
# plotting
#

style.use('dark_background')

# isotherms tilted
fig = plt.figure(figsize=(7,7.))
ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
ax.spines['right'].set_color('none')
#ax.spines['top'].set_color('none')
ax.spines['left'].set_color('none')
ax.plot(phi1,z,'--',linewidth=lw1+2,alpha=aph,color=color1)
ax.plot(phi2,z,linewidth=lw1+2,alpha=aph,color=color1)

ax.text(.09,.02,'z = 0',fontsize=20.)
ax.text(.09,1.02,'z = H',fontsize=20.)
plt.yticks([])

plt.xlabel('Amplitude')

plt.savefig('figs/sqg_waves',bbox_inches='tight',transparent=True)    



