
import numpy as np
import matplotlib.pyplot as plt
from mpltools import style


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

    fig = plt.figure(facecolor='black',figsize=(12,7.))
    ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    ax.plot([.12,.72],[.3,.3],'--',color='w',linewidth=.1,alpha=.65)

    ax.contour(yi,zi,temp,np.linspace(-.15,0.,4),cmap='Spectral_r',
            linewidths=5.)
    #plt.contourf(yi,zi,temp,np.linspace(-.15,0.,6),cmap='Spectral_r')
    #ax.text(0.05,.8,'Buoyant',fontsize=20.)
    #ax.arrow( 0.2, 0., -.1, 0.0, fc="k", ec="k",head_width=0.03
    #        , head_length=0.05 )

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

    ax.arrow( 0.8, 0.15, .1, 0.0, color='w',fc="w", ec="w",head_width=0.05
          , head_length=0.034 )

    # water parcels
    ax.plot([.4,.72],[.3,.6],color='w',alpha=.85)

    ax.plot(.4,.3,'o',markersize=14.,color='#ffd700')
    ax.plot(.4,.7,'o',markersize=14.,color='#b22222')
    ax.plot(.72,.6,'o',markersize=14.,color='#9acd32')


    plt.axis('off')

ax.text(.17,.34,u'$\phi$',fontsize=24,color='w')
ax.text(.47,.335,u'$\gamma$',fontsize=24,color='w')

ann = ax.annotate("",
              xy=(0.15364583333331, 0.301020408163265), xycoords='data',
              xytext=(0.143229166666, 0.3622448979591), textcoords='data',
              size=20, va="center", ha="center",
              bbox=dict(boxstyle="round4", fc="w"),
              arrowprops=dict(arrowstyle="-",
                              connectionstyle="arc3,rad=-.6",
                              relpos=(0., 0.),
                              fc="w"),fontsize=25. 
              )

ann = ax.annotate("",
              xy=(0.4472916666666, 0.2959183673469), xycoords='data',
              xytext=(0.44078125, 0.35204081632), textcoords='data',
              size=20, va="center", ha="center",
              bbox=dict(boxstyle="round4", fc="w"),
              arrowprops=dict(arrowstyle="-",
                              connectionstyle="arc3,rad=-.6",
                              relpos=(0., 0.),
                              fc="w",),fontsize=25.
              )

ax.text(.4,.23,'A',fontsize=20.,color='w')
ax.text(.4,.63,'B',fontsize=20.,color='w')
ax.text(.72,.53,'C',fontsize=20.,color='w')


ax.text(0.0,0.82,u'$f/2$',fontsize=23.,color='w')
ax.text(0.05,0.92,u'$z$',fontsize=25.,color='w')

ax.text(0.8,.05,'Poleward',fontsize=20.,color='w')
    
plt.savefig('figs/parcels',bbox_inches='tight',transparent=True)    

    
# parcels 2: change A with B
with plt.xkcd():

    fig = plt.figure(facecolor='black',figsize=(12,7.))
    ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    ax.plot([.12,.72],[.3,.3],'--',color='w',linewidth=.1,alpha=.65)

    ax.contour(yi,zi,temp,np.linspace(-.15,0.,4),cmap='Spectral_r',
            linewidths=5.)
    #plt.contourf(yi,zi,temp,np.linspace(-.15,0.,6),cmap='Spectral_r')
    #ax.text(0.05,.8,'Buoyant',fontsize=20.)
    #ax.arrow( 0.2, 0., -.1, 0.0, fc="k", ec="k",head_width=0.03
    #        , head_length=0.05 )

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

    ax.arrow( 0.8, 0.15, .1, 0.0, color='w',fc="w", ec="w",head_width=0.05
          , head_length=0.034 )

    # water parcels
    ax.plot([.4,.72],[.3,.6],color='w',alpha=.85)

    ax.plot(.4,.7,'o',markersize=14.,color='#ffd700')
    ax.plot(.4,.3,'o',markersize=14.,color='#b22222')
    ax.plot(.72,.6,'o',markersize=14.,color='#9acd32')


    plt.axis('off')

ax.text(.17,.34,u'$\phi$',fontsize=24,color='w')
ax.text(.47,.335,u'$\gamma$',fontsize=24,color='w')

ann = ax.annotate("",
              xy=(0.15364583333331, 0.301020408163265), xycoords='data',
              xytext=(0.143229166666, 0.3622448979591), textcoords='data',
              size=20, va="center", ha="center",
              bbox=dict(boxstyle="round4", fc="w"),
              arrowprops=dict(arrowstyle="-",
                              connectionstyle="arc3,rad=-.6",
                              relpos=(0., 0.),
                              fc="w"),fontsize=25. 
              )

ann = ax.annotate("",
              xy=(0.4472916666666, 0.2959183673469), xycoords='data',
              xytext=(0.44078125, 0.35204081632), textcoords='data',
              size=20, va="center", ha="center",
              bbox=dict(boxstyle="round4", fc="w"),
              arrowprops=dict(arrowstyle="-",
                              connectionstyle="arc3,rad=-.6",
                              relpos=(0., 0.),
                              fc="w",),fontsize=25.
              )

ax.text(.4,.23,'B',fontsize=20.,color='w')
ax.text(.4,.63,'A',fontsize=20.,color='w')
ax.text(.72,.53,'C',fontsize=20.,color='w')


ax.text(0.0,0.82,u'$f/2$',fontsize=23.,color='w')
ax.text(0.05,0.92,u'$z$',fontsize=25.,color='w')

ax.text(0.8,.05,'Poleward',fontsize=20.,color='w')
    
plt.savefig('figs/parcels2',bbox_inches='tight',transparent=True)    

# now change A with C 
with plt.xkcd():

    fig = plt.figure(facecolor='black',figsize=(12,7.))
    ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    ax.plot([.12,.72],[.3,.3],'--',color='w',linewidth=.1,alpha=.65)

    ax.contour(yi,zi,temp,np.linspace(-.15,0.,4),cmap='Spectral_r',
            linewidths=5.)
    #plt.contourf(yi,zi,temp,np.linspace(-.15,0.,6),cmap='Spectral_r')
    #ax.text(0.05,.8,'Buoyant',fontsize=20.)
    #ax.arrow( 0.2, 0., -.1, 0.0, fc="k", ec="k",head_width=0.03
    #        , head_length=0.05 )

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

    ax.arrow( 0.8, 0.15, .1, 0.0, color='w',fc="w", ec="w",head_width=0.05
          , head_length=0.034 )

    # water parcels
    ax.plot([.4,.72],[.3,.6],color='w',alpha=.85)

    ax.plot(.72,.6,'o',markersize=14.,color='#ffd700')
    ax.plot(.4,.7,'o',markersize=14.,color='#b22222')
    ax.plot(.4,.3,'o',markersize=14.,color='#9acd32')


    plt.axis('off')

ax.text(.17,.34,u'$\phi$',fontsize=24,color='w')
ax.text(.47,.335,u'$\gamma$',fontsize=24,color='w')

ann = ax.annotate("",
              xy=(0.15364583333331, 0.301020408163265), xycoords='data',
              xytext=(0.143229166666, 0.3622448979591), textcoords='data',
              size=20, va="center", ha="center",
              bbox=dict(boxstyle="round4", fc="w"),
              arrowprops=dict(arrowstyle="-",
                              connectionstyle="arc3,rad=-.6",
                              relpos=(0., 0.),
                              fc="w"),fontsize=25. 
              )

ann = ax.annotate("",
              xy=(0.4472916666666, 0.2959183673469), xycoords='data',
              xytext=(0.44078125, 0.35204081632), textcoords='data',
              size=20, va="center", ha="center",
              bbox=dict(boxstyle="round4", fc="w"),
              arrowprops=dict(arrowstyle="-",
                              connectionstyle="arc3,rad=-.6",
                              relpos=(0., 0.),
                              fc="w",),fontsize=25.
              )

ax.text(.4,.23,'C',fontsize=20.,color='w')
ax.text(.4,.63,'B',fontsize=20.,color='w')
ax.text(.72,.53,'A',fontsize=20.,color='w')


ax.text(0.0,0.82,u'$f/2$',fontsize=23.,color='w')
ax.text(0.05,0.92,u'$z$',fontsize=25.,color='w')

ax.text(0.8,.05,'Poleward',fontsize=20.,color='w')
    
plt.savefig('figs/parcels3',bbox_inches='tight',transparent=True)  
