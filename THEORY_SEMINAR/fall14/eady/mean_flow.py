
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection


plt.close('all')
plt.rcParams.update({'font.size': 30, 'legend.handlelength'  : 1.5
        , 'legend.markerscale': 14., 'legend.linewidth': 3.})


y = np.linspace(0.,1.,30)
z = np.linspace(0.,1.,30)

yi,zi = np.meshgrid(y,z)

A,B = .25,.1  

temp = -A*yi + B*zi


#
# plotting
#

#plt.style.use('dark_background')


# isotherms tilted

with plt.xkcd():
    fig = plt.figure(facecolor='black',figsize=(12,9.))
    ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    ax.plot([-.1,1.2],[0.,0.],color='w',linewidth=3.,alpha=.65)
    ax.plot([-.1,1.2],[1.,1.],color='w',linewidth=3.,alpha=.65)

    ax.contour(yi,zi,temp,np.linspace(-.15,0.,4),cmap='Spectral_r',
            linewidths=5.)

    ax.arrow( 0.1, 0.7, 0, 0.2, fc="w", ec="w",head_width=0.03
            , head_length=0.05,)

    ann = ax.annotate("",
                  xy=(0.15, 0.86), xycoords='data',
                  xytext=(0.05, 0.85), textcoords='data',
                  size=20, va="center", ha="center",
                  bbox=dict(boxstyle="round4", fc="w"),
                  arrowprops=dict(arrowstyle="-|>",
                                  connectionstyle="arc3,rad=.6",
                                  relpos=(0., 0.),
                                  fc="w"), 
                  )
    ax.text(0.0,0.82,u'f/2',fontsize=20.)

    ax.text(0.8,.05,'Poleward',fontsize=20.,color='k')
    ax.arrow( 0.8, 0.15, .1, 0.0, fc="w", ec="w",head_width=0.05
            , head_length=0.034 )

    plt.axis('off')


    # shear
    plt.savefig('figs/mean_flow',bbox_inches='tight',transparent=True)    

    fig = plt.figure(facecolor='black',figsize=(12,9.))
    ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    ax.plot([-.1,1.2],[0.,0.],color='w',linewidth=3.,alpha=.65)
    ax.plot([-.1,1.2],[1.,1.],color='w',linewidth=3.,alpha=.65)

    ax.contour(yi,zi,temp,np.linspace(-.15,0.,4),cmap='Spectral_r',
            linewidths=5.)

    ax.arrow( 0.1, 0.7, 0, 0.2, fc="w", ec="w",head_width=0.03
            , head_length=0.05,)

    ann = ax.annotate("",
                  xy=(0.15, 0.86), xycoords='data',
                  xytext=(0.05, 0.85), textcoords='data',
                  size=20, va="center", ha="center",
                  bbox=dict(boxstyle="round4", fc="w"),
                  arrowprops=dict(arrowstyle="-|>",
                                  connectionstyle="arc3,rad=.6",
                                  relpos=(0., 0.),
                                  fc="w"), 
                  )
    ax.text(0.0,0.82,u'f/2',fontsize=20.)

    ax.text(0.8,.05,'Poleward',fontsize=20.,color='k')
    ax.arrow( 0.8, 0.15, .1, 0.0, fc="w", ec="w",head_width=0.05
            , head_length=0.034 )


    ax.text(-.09,.02,'z = 0',fontsize=20.)
    ax.text(-.09,1.02,'z = H',fontsize=20.)

    # constant shear
    circle1=plt.Circle((.6,.8),.14,edgecolor='w',facecolor='k')
    fig.gca().add_artist(circle1) 
    ax.plot(.6,.8,'o',color='w',markersize=10.)

    circle2=plt.Circle((.6,.53),.11,edgecolor='w',facecolor='k')
    fig.gca().add_artist(circle2) 
    ax.plot(.6,.53,'o',color='w',markersize=10.)

    circle3=plt.Circle((.6,.33),.07,edgecolor='w',facecolor='k')
    fig.gca().add_artist(circle3) 
    ax.plot(.6,.33,'o',color='w',markersize=8.)

    circle4=plt.Circle((.6,.19),.05,edgecolor='w',facecolor='k')
    fig.gca().add_artist(circle4) 
    ax.plot(.6,.19,'o',color='w',markersize=7.)

    circle5=plt.Circle((.6,.097),.025,edgecolor='w',facecolor='k')
    fig.gca().add_artist(circle5) 
    ax.plot(.6,.097,'o',color='w',markersize=5.)

    circle6=plt.Circle((.6,.035),.01,edgecolor='w',facecolor='k')
    fig.gca().add_artist(circle6) 
    ax.plot(.6,.035,'o',color='w',markersize=1.)

    plt.axis('off')

    plt.savefig('figs/mean_flow_shear',bbox_inches='tight',transparent=True)    




