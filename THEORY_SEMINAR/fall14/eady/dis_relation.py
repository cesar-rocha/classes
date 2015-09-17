
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from mpltools import style


plt.close('all')
plt.rcParams.update({'font.size': 30, 'legend.handlelength'  : 1.5
        , 'legend.markerscale': 14., 'legend.linewidth': 3.})


mu2 = np.linspace(0.,2.,)

tmu2 = sp.tanh(mu2)
cmu2 = 1./tmu2

color1 = '#ffd700'
color3 = '#9acd32'
color2 = '#ff0000'

lw1=4
aph=.7

# set the linewidth of each legend object
def leg_width(lg,fs):
    for legobj in lg.legendHandles:
        legobj.set_linewidth(fs)

#
# plotting
#

style.use('dark_background')



# isotherms tilted

fig = plt.figure(figsize=(17,9.))
ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')


ax.plot(mu2,mu2,linewidth=lw1,alpha=aph,color=color2)
ax.plot(mu2,cmu2,linewidth=lw1,alpha=aph,color=color1)
ax.plot(mu2,tmu2,linewidth=lw1,alpha=aph,color=color3)

ax.text(.35,1.6,u'$coth(\mu/2)$',fontsize=25.,color=color1)
ax.text(1.45,.8,u'$tanh(\mu/2)$',fontsize=25.,color=color3)
ax.text(1.7,1.6,u'$\mu/2$',fontsize=25.,color=color2)

plt.ylim((0.,2.))
plt.yticks([])

ax.set_xticks([.5,1.2])
ax.set_xticklabels([r'$\mu$ = 1 ', r'$\mu_c$ $\approx$ 2.4'])
#plt.axis('off')
plt.savefig('figs/tanh_coth',bbox_inches='tight',transparent=True)    

fig = plt.figure(figsize=(17,9.))
ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')


ax.plot([1.2,1.2],[0.,5.],'--',color='w',linewidth=2.)

ax.plot(mu2,mu2,linewidth=lw1,alpha=aph,color=color2)
ax.plot(mu2,cmu2,linewidth=lw1,alpha=aph,color=color1)
ax.plot(mu2,tmu2,linewidth=lw1,alpha=aph,color=color3)

ax.text(.35,1.6,u'$coth(\mu/2)$',fontsize=32.,color=color1)
ax.text(1.45,.75,u'$tanh(\mu/2)$',fontsize=32.,color=color3)
ax.text(1.7,1.55,u'$\mu/2$',fontsize=32.,color=color2)


ax.text(0.45, 1.2, "Unstable", size=25, rotation=0.,
         ha="center", va="center",
         bbox = dict(boxstyle="round",ec='w',fc=color2))

ax.text(1.65,.35, "Stable", size=25, rotation=0.,
         ha="center", va="center",
         bbox = dict(boxstyle="round",ec='w',fc='b'))

plt.ylim((0.,2.))
plt.yticks([])

ax.set_xticks([.5,1.2])
ax.set_xticklabels([r'$\mu$ = 1 ', r'$\mu_c$ $\approx$ 2.4'])
#plt.axis('off')
plt.savefig('figs/tanh_coth_unstable',bbox_inches='tight',transparent=True)    


#
# growth rate
#

k = np.linspace(0.01,3.01,100.)

def grate(k,L):
    l2 = (np.pi**2)/(4*(L**2))
    mu = np.sqrt(k**2 + l2) 
    th = sp.tanh(mu/2.)
    cth = 1./th
    sigma = (k/mu)*np.sqrt( (mu/2.-cth)*(th-mu/2.)  )
    sigma[np.isnan(sigma)] = 0.
    return sigma

sigma = grate(k,1.)
sigma2 = grate(k,2.)
sigma4 = grate(k,4.)
sigma8 = grate(k,8.)

fig = plt.figure(figsize=(17,9.))
ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')


ax.plot(k,sigma,linewidth=5.,alpha=.7,label='1')
ax.plot(k,sigma2,linewidth=5.,alpha=.7,label='2')
ax.plot(k,sigma4,linewidth=5.,alpha=.7,label='4')
ax.plot(k,sigma8,linewidth=5.,alpha=.7,label='8')

plt.xlabel(r'$k L_d$',fontsize=30.)
plt.ylabel(r'$\sigma L_d/\Lambda H$',fontsize=30.)

lg = plt.legend(loc=1,title= r'$L/L_d$',fontsize=30., prop={'size':35}, numpoints=1)
lg.draw_frame(False)
leg_width(lg,8.)

plt.xlim(0.,3.)
plt.savefig('figs/sigma',bbox_inches='tight',transparent=True)    

#
# ci, cr
#
k = np.linspace(0.01,5.01,1000.)
L = 8.
l2 = (np.pi**2)/(4*(L**2))
mu = np.sqrt(k**2 + l2)
th = sp.tanh(mu/2.)
cth = 1./th

cr =  np.zeros(k.size)
cr2 = np.zeros(k.size)

ci = np.copy(cr)
ci2 = np.copy(cr)

ki,ks=(k<=2.4),(k>=2.4)

cr[ki],cr2[ki] = .5, .5
cr[ks] = .5 + (1./mu[ks])*np.sqrt( (mu[ks]/2. - cth[ks])*(mu[ks]/2. - th[ks]) )
cr2[ks] = .5 - (1./mu[ks])*np.sqrt( (mu[ks]/2. - cth[ks])*(mu[ks]/2. - th[ks]) )

ci[ki] =   (1./mu[ki])*np.sqrt( ( cth[ki] - mu[ki]/2.)*(mu[ki]/2. - th[ki]) )
ci2[ki]= - (1./mu[ki])*np.sqrt( ( cth[ki] - mu[ki]/2.)*(mu[ki]/2. - th[ki]) )


fig = plt.figure(figsize=(17,9.))

ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

ax.plot([2.4,2.4],[-.4,.9],'--',color='w')
ax.text(2.395,.92,r'$\mu_c$',color='w')

ax.plot(k,0.*k,'--',color='w')
ax.plot(k,ci,linewidth=3.,color=color2,alpha=.7,label='c_i')
ax.plot(k,ci2,linewidth=3.,color=color2,alpha=.7)
ax.plot(k,cr,linewidth=3.,color=color1,alpha=.7,label='c_r')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.plot(k,cr2,linewidth=3.,color=color1,alpha=.7)

plt.xlabel(r'$\mu$')
plt.ylabel(r'Wave speed/$\Lambda H$')

lg = plt.legend(loc=2,title= r'', prop={'size':22}, numpoints=1)
lg.draw_frame(False)
leg_width(lg,8.)

plt.xlim(0.,3.)
plt.savefig('figs/speeds',bbox_inches='tight',transparent=True)    


#
# vertical structure of most unstable wave
#

z = np.linspace(0.,1.,100)

sigma8 = grate(k,8.)
ikmax = sigma8.argmax()
mu = mu[ikmax]
co = sp.cosh(mu*z)
si = sp.sinh(mu*z)

cr,ci = .5,ci[ikmax]
c2 = cr**2 + ci**2

co = sp.cosh(mu*z)
si = sp.sinh(mu*z)
phi_r = co - (cr/(mu*c2))*si
phi_i = (ci/(mu*c2))*si

phi_abs = np.sqrt( phi_r**2 + phi_i**2 )

phi_phase = np.zeros(z.size)
for i in range(z.size):
    phi_phase[i] = sp.math.atan2( phi_i[i],phi_r[i] )

fig = plt.figure(figsize=(12,12.))
ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
ax.plot(phi_abs,z,linewidth=7.)
plt.xlabel('Amplitude')
plt.ylabel('z/H')
plt.xlim(0.5,1.)
plt.savefig('figs/amplitude',bbox_inches='tight',transparent=True)    

fig = plt.figure(figsize=(12,12.))
ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
ax.plot(phi_phase,z,linewidth=7.)
plt.xlabel('Phase')
plt.ylabel('z/H')

ax.set_xticks([0.,np.pi/4.,np.pi/2.])
ax.set_xticklabels(['0',r'$\pi/4$',r'$\pi/2$'])

plt.xlim(0.,np.pi/2.)
plt.savefig('figs/phase',bbox_inches='tight',transparent=True)    

#
# wave structure
#

k=k[ikmax]

x = np.linspace(0.,3.*np.pi,100)

xi,zi = np.meshgrid(x,z)

dz = 1./z.size
phi_abs_z = np.gradient(phi_abs,dz) 
phi_phase_z = np.gradient(phi_phase,dz) 

phi_abs = np.repeat(phi_abs,100).reshape(100,100)
phi_phase = np.repeat(phi_phase,100).reshape(100,100)

phi_abs_z = np.repeat(phi_abs_z,100).reshape(100,100)
phi_phase_z = np.repeat(phi_phase_z,100).reshape(100,100)

# solutions
t = 0.
cr = .5
R = 1.
phase_0 = -np.pi/4

psi = R*phi_abs*np.cos(x - cr*t + phi_phase -phase_0)
v   = -R*phi_abs*np.sin(x - cr*t + phi_phase -phase_0)
b = R*( phi_abs_z*np.cos(x - cr*t + phi_phase -phase_0) - 
        phi_abs*phi_phase_z*np.sin(x - cr*t + phi_phase -phase_0))  
w = -R*( cr*phi_abs_z*np.sin(x - cr*t + phi_phase -phase_0) + 
        cr*phi_abs*phi_phase_z*np.cos(x - cr*t + phi_phase -phase_0))  


fig = plt.figure(figsize=(15,10.))

ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
ax.contour(xi,zi,psi,10,colors='w',linewidths=3.)
ax.set_xticks([0.,np.pi,2.*np.pi,3.*np.pi])
ax.set_xticklabels(['0',r'$\pi$',r'$2\pi$',r'$3\pi$'])
plt.xlim(0.,3.*np.pi)

plt.xlabel(r'$x k$')
plt.ylabel(r'$z/H$')
plt.title(r'Stream function',fontsize=30.)

plt.savefig('figs/psi',bbox_inches='tight',transparent=True)    

fig = plt.figure(figsize=(15,10.))
ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
ax.contour(xi,zi,v,10,colors='w',linewidths=3.)
ax.set_xticks([0.,np.pi,2.*np.pi,3.*np.pi])
ax.set_xticklabels(['0',r'$\pi$',r'$2\pi$',r'$3\pi$'])
plt.xlim(0.,3.*np.pi)
plt.xlabel(r'$x k$')
plt.ylabel(r'$z/H$')
plt.title(r'Meridional velocity',fontsize=30.)
plt.savefig('figs/v',bbox_inches='tight',transparent=True) 

fig = plt.figure(figsize=(15,10.))
ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
ax.contour(xi,zi,w,10,colors='w',linewidths=3.)
ax.set_xticks([0.,np.pi,2.*np.pi,3.*np.pi])
ax.set_xticklabels(['0',r'$\pi$',r'$2\pi$',r'$3\pi$'])
plt.xlim(0.,3.*np.pi)
plt.xlabel(r'$x k$')
plt.ylabel(r'$z/H$')
plt.title(r'Vertical velocity',fontsize=30.)
plt.savefig('figs/w',bbox_inches='tight',transparent=True) 


blue = '#6495ed'
red = '#ff6347'

fig = plt.figure(figsize=(15,10.))
ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))

ax.contour(xi,zi,b,np.linspace(0.1,b.max(),5),colors=red,linewidths=3.)
ax.contour(xi,zi,-b,np.linspace(0.1,np.abs(b.min()),5),colors=blue,linewidths=3.)

ax.set_xticks([0.,np.pi,2.*np.pi,3.*np.pi])
ax.set_xticklabels(['0',r'$\pi$',r'$2\pi$',r'$3\pi$'])
plt.xlim(0.,3.*np.pi)
plt.xlabel(r'$x k$')
plt.ylabel(r'$z/H$')
plt.title(r'Buoyancy',fontsize=30.)
plt.savefig('figs/b',bbox_inches='tight',transparent=True) 

#
# correlations
#

fig = plt.figure(figsize=(15,10.))
ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))

ax.contour(xi,zi,b,np.linspace(0.1,b.max(),5),colors=red,linewidths=3.)
ax.contour(xi,zi,-b,np.linspace(0.1,np.abs(b.min()),5),colors=blue,linewidths=3.)
ax.contour(xi,zi,v,10,colors='w',linewidths=3.)
ax.set_xticks([0.,np.pi,2.*np.pi,3.*np.pi])
ax.set_xticklabels(['0',r'$\pi$',r'$2\pi$',r'$3\pi$'])
plt.xlim(0.,3.*np.pi)
plt.xlabel(r'$x k$')
plt.ylabel(r'$z/H$')
plt.title(r'$\overline{vb}>0$',fontsize=30.)
plt.savefig('figs/bv',bbox_inches='tight',transparent=True) 

fig = plt.figure(figsize=(15,10.))
ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))

ax.contour(xi,zi,b,np.linspace(0.1,b.max(),5),colors=red,linewidths=3.)
ax.contour(xi,zi,-b,np.linspace(0.1,np.abs(b.min()),5),colors=blue,linewidths=3.)
ax.contour(xi,zi,w,10,colors='w',linewidths=3.)
ax.set_xticks([0.,np.pi,2.*np.pi,3.*np.pi])
ax.set_xticklabels(['0',r'$\pi$',r'$2\pi$',r'$3\pi$'])
plt.xlim(0.,3.*np.pi)
plt.xlabel(r'$x k$')
plt.ylabel(r'$z/H$')
plt.title(r'$\overline{wb}>0$',fontsize=30.)
plt.savefig('figs/bw',bbox_inches='tight',transparent=True) 


#
# vertical structure stable modes
#

mu = 5.
phi_upper = np.exp(mu*(zi-1))
phi_lower = np.exp(-mu*zi)

phi_phase = 0.

phi_upper_z = mu*phi_upper
phi_lower_z = -mu*phi_lower

psi = R*(phi_upper*np.cos(x - cr*t + phi_phase -phase_0) +  phi_lower*np.cos(x - cr*t + phi_phase -phase_0)) 

v = -R*(phi_upper*np.sin(x - cr*t + phi_phase -phase_0) +  phi_lower*np.sin(x - cr*t + phi_phase -phase_0)) 

b =  -R*(phi_upper_z*phi_upper*np.sin(x - cr*t + phi_phase -phase_0) +  phi_lower_z*phi_lower*np.sin(x - cr*t + phi_phase -phase_0)) 




fig = plt.figure(figsize=(15,10.))

ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
ax.contour(xi,zi,psi,10,colors='w',linewidths=3.)
ax.set_xticks([0.,np.pi,2.*np.pi,3.*np.pi])
ax.set_xticklabels(['0',r'$\pi$',r'$2\pi$',r'$3\pi$'])
plt.xlim(0.,3.*np.pi)

plt.xlabel(r'$x k$')
plt.ylabel(r'$z/H$')
plt.title(r'Stream function',fontsize=30.)

plt.savefig('figs/psi_stable',bbox_inches='tight',transparent=True)    

fig = plt.figure(figsize=(15,10.))
ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
ax.contour(xi,zi,v,10,colors='w',linewidths=3.)
ax.set_xticks([0.,np.pi,2.*np.pi,3.*np.pi])
ax.set_xticklabels(['0',r'$\pi$',r'$2\pi$',r'$3\pi$'])
plt.xlim(0.,3.*np.pi)
plt.xlabel(r'$x k$')
plt.ylabel(r'$z/H$')
plt.title(r'Meridional velocity',fontsize=30.)
plt.savefig('figs/v_stable',bbox_inches='tight',transparent=True) 

#fig = plt.figure(figsize=(15,10.))
#ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
#ax.contour(xi,zi,w,10,colors='w',linewidths=3.)
#ax.set_xticks([0.,np.pi,2.*np.pi,3.*np.pi])
#ax.set_xticklabels(['0',r'$\pi$',r'$2\pi$',r'$3\pi$'])
#plt.xlim(0.,3.*np.pi)
#plt.xlabel(r'$x k$')
#plt.ylabel(r'$z/H$')
#plt.title(r'Vertical velocity',fontsize=30.)
#plt.savefig('figs/w',bbox_inches='tight',transparent=True) 

fig = plt.figure(figsize=(15,10.))
ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))

ax.contour(xi,zi,b,np.linspace(0.5,b.max(),5),colors=red,linewidths=3.)
ax.contour(xi,zi,-b,np.linspace(0.5,np.abs(b.min()),5),colors=blue,linewidths=3.)

ax.set_xticks([0.,np.pi,2.*np.pi,3.*np.pi])
ax.set_xticklabels(['0',r'$\pi$',r'$2\pi$',r'$3\pi$'])
plt.xlim(0.,3.*np.pi)
plt.xlabel(r'$x k$')
plt.ylabel(r'$z/H$')
plt.title(r'Buoyancy',fontsize=30.)
plt.savefig('figs/b_stable',bbox_inches='tight',transparent=True) 


#
# correlations
#


fig = plt.figure(figsize=(15,10.))
ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))

ax.contour(xi,zi,b,np.linspace(0.5,b.max(),5),colors=red,linewidths=3.)
ax.contour(xi,zi,-b,np.linspace(0.5,np.abs(b.min()),5),colors=blue,linewidths=3.)
ax.contour(xi,zi,v,10,colors='w',linewidths=3.)
ax.set_xticks([0.,np.pi,2.*np.pi,3.*np.pi])
ax.set_xticklabels(['0',r'$\pi$',r'$2\pi$',r'$3\pi$'])
plt.xlim(0.,3.*np.pi)
plt.xlabel(r'$x k$')
plt.ylabel(r'$z/H$')
plt.title(r'$\overline{vb}=0$',fontsize=30.)
plt.savefig('figs/bv_stable',bbox_inches='tight',transparent=True) 

#fig = plt.figure(figsize=(15,10.))
#ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
#
#ax.contour(xi,zi,b,np.linspace(0.1,b.max(),5),colors=red,linewidths=3.)
#ax.contour(xi,zi,-b,np.linspace(0.1,np.abs(b.min()),5),colors=blue,linewidths=3.)
#ax.contour(xi,zi,w,10,colors='w',linewidths=3.)
#ax.set_xticks([0.,np.pi,2.*np.pi,3.*np.pi])
#ax.set_xticklabels(['0',r'$\pi$',r'$2\pi$',r'$3\pi$'])
#plt.xlim(0.,3.*np.pi)
#plt.xlabel(r'$x k$')
#plt.ylabel(r'$z/H$')
#plt.title(r'$\overline{wb}=0$',fontsize=30.)
#plt.savefig('figs/bw',bbox_inches='tight',transparent=True) 


#
# The sear
#

z = np.linspace(0.,1.,12)

U = z
W = 0
x = z*0

fig = plt.figure(figsize=(8,12.))
ax = fig.add_axes((0.1, 0.2, 0.8, 0.7))
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.spines['left'].set_color('none')
ax.quiver(x,z,10*U, W,color='w',linewidths=(1,), edgecolors=('k'), headaxislength=5,scale=15)

plt.xlim(0.,.25)
plt.ylim(0.,1.1)
plt.xticks([])
plt.yticks([])
plt.savefig('figs/shear_vecs',bbox_inches='tight',transparent=True)
