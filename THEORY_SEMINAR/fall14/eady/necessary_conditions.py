
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.axislines import SubplotZero
import matplotlib.lines as lines

from mpltools import style


plt.close('all')
plt.rcParams.update({'font.size': 30, 'legend.handlelength'  : 1.5
        , 'legend.markerscale': 14., 'legend.linewidth': 3.})


color1 = '#ffd700'
color3 = '#9acd32'
color2 = '#ff0000'

lw1=4
aph=.7

style.use('dark_background')

#
# pv gradient reversal
#

z = np.linspace(0.,np.pi,30)
qy = -np.cos(z)

def make_yaxis(ax, xloc=0, offset=0.05, **props):
    ymin, ymax = ax.get_ylim()
    locs = [loc for loc in ax.yaxis.get_majorticklocs()
            if loc>=ymin and loc<=ymax]
    tickline, = ax.plot([xloc]*len(locs), locs, linestyle='',
            marker=lines.TICKLEFT, **props)
    axline, = ax.plot([xloc, xloc], [ymin, ymax], **props)
    tickline.set_clip_on(False)
    axline.set_clip_on(False)

    for loc in locs:
        ax.text(xloc-offset, loc, '%1.1f'%loc,
                verticalalignment='center',
                horizontalalignment='right')



with plt.xkcd():
    fig = plt.figure(figsize=(8,10.))
    ax = fig.add_subplot(111)
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    plt.ylim(0,np.pi)
    ax.set_yticks([])
    ax.set_xticks([0])
    ax.set_xticklabels(['0'])
    plt.xlim(-1.6,1.6)
    c1 = '#6495ed'
    c2 = color='#ff6347'

    make_yaxis(ax, 0, offset=10,linewidth=2)

ax.fill_betweenx(z[0:15],qy[0:15],0,color=c1,alpha=.3)
ax.fill_betweenx(z[15:],qy[15:],0,color=c2,alpha=.3)
ax.text(-1.57,0.05,r'z = 0')
ax.text(-1.57,np.pi+0.05,r'z = H')

plt.text(-.65,np.pi/4-.2,r'$Q_y < 0$',color='w')
plt.text(.18,3*np.pi/4+.1,r'$Q_y > 0$',color='w')
plt.savefig('figs/pv_reversal',bbox_inches='tight',transparent=True)    


qy = np.ones(z.size)
with plt.xkcd():
    fig = plt.figure(figsize=(8,10.))
    ax = fig.add_subplot(111)
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    plt.ylim(0,np.pi)
    ax.set_yticks([])
    ax.set_xticks([0])
    ax.set_xticklabels(['0'])
    plt.xlim(-1.6,1.6)
    c1 = '#6495ed'
    c2 = color='#ff6347'

    ax.plot([-1.7,1.7],[np.pi,np.pi],color='b',linewidth=40)

    make_yaxis(ax, 0, offset=10,linewidth=2)

ax.fill_betweenx(z,qy,0,color=c2,alpha=.3)
ax.text(-1.57,0.05,r'z = 0')
ax.text(-1.57,np.pi+0.05,r'z = H')
ax.text(.15,np.pi/2.,r'$Q_y > 0$',color='w')
ax.text(1.,np.pi+0.05,r'$U_z < 0$',color='w')
plt.savefig('figs/uz_upper',bbox_inches='tight',transparent=True)    

\int_0^L dy \bigg\{ \int_0^H \frac{|\phi^2|}{|U-c|^2} Q_y dz + \frac{f_0^2}{N^2(z)} \frac{|\phi^2|}{|U-c|^2} U_z  \bigg|_0^H\bigg\} = 0


qy = np.ones(z.size)
with plt.xkcd():
    fig = plt.figure(figsize=(8,10.))
    ax = fig.add_subplot(111)
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    plt.ylim(0,np.pi)
    ax.set_yticks([])
    ax.set_xticks([0])
    ax.set_xticklabels(['0'])
    plt.xlim(-1.6,1.6)
    c1 = '#6495ed'
    c2 = color='#ff6347'

    ax.plot([-1.7,1.7],[0.,0.],color='r',linewidth=40)

    make_yaxis(ax, 0, offset=10,linewidth=2)
    ax.plot([-1.7,1.7],[0,0],color='r',linewidth=40)

ax.fill_betweenx(z,qy,0,color=c2,alpha=.3)

ax.text(-1.57,0.05,r'z = 0')
ax.text(-1.57,np.pi+0.05,r'z = H')
ax.text(.15,np.pi/2.,r'$Q_y > 0$',color='w')
ax.text(1.,0.05,r'$U_z > 0$',color='w')
plt.savefig('figs/uz_lower',bbox_inches='tight',transparent=True)    

with plt.xkcd():
    fig = plt.figure(figsize=(8,10.))
    ax = fig.add_subplot(111)
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    plt.ylim(0,np.pi)
    ax.set_yticks([])
    ax.set_xticks([0])
    ax.set_xticklabels(['0'])
    plt.xlim(-1.6,1.6)
    c1 = '#6495ed'
    c2 = color='#ff6347'


    make_yaxis(ax, 0, offset=10,linewidth=2)

    ax.plot([-1.7,1.7],[np.pi,np.pi],color='r',linewidth=40)
    ax.plot([-1.7,1.7],[0,0],color='r',linewidth=40)

ax.text(-1.57,0.05,r'z = 0')
ax.text(-1.57,np.pi+0.05,r'z = H')
ax.text(1.,np.pi+0.05,r'$U_z > 0$',color='w')
ax.text(1.,0.05,r'$U_z > 0$',color='w')
plt.savefig('figs/uz_both',bbox_inches='tight',transparent=True)    

