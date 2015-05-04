# MAE290C, Homework 2, Problem 3
#   Plotting script
# Cesar B. Rocha

import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import sys

plt.close('all')

N = int(sys.argv[1])

if N<100:
    Ns = "00"+str(N)
elif N < 1000:
    Ns = "0"+str(N)
else:
    Ns = str(N)


time = np.loadtxt('output/Burgers'+Ns+'.time')
ini = np.loadtxt('output/Burgers'+Ns+'.ini')
phys = np.loadtxt('output/Burgers'+Ns+'.phys')
four = np.loadtxt('output/Burgers'+Ns+'.four')
k = np.loadtxt('output/Burgers'+Ns+'.k')

ix,jx = phys.shape

tplot = np.array([0.25,.5,1.,2.,5.])
iplot = np.empty_like(tplot)
lb = np.empty_like(tplot)

for i in range(tplot.size):
    iplot[i] = ((np.abs(tplot[i]-time)).argmin())

# plot soln in physical odmain
fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(111)
ax.plot(ini[:,0],ini[:,1], label='0')
for i in range(iplot.size):
    ax.plot(ini[:,0],phys[iplot[i],:], label=str(tplot[i]))

lg = plt.legend(loc=1)
ax.text(.1,1.4,'N = '+ str(N),fontsize=16)

ax.set_xticks(np.array([0.,pi/2,pi,3*pi/2,2*pi]))
ax.set_xticklabels(np.array([r'$0$',r'$\pi/2$',r'$\pi$',r'$3\pi/2$',r'$2\pi$']))
ax.set_xlim(0,2*pi)
ax.set_ylim(-1.6,1.6)

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$u$')

plt.savefig('figs/burgers_phys_'+str(N))

# plot spectrum
fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(111)
ax.semilogx(k,np.abs(np.fft.rfft(ini[:,1]))/N, label='0')
for i in range(iplot.size):
    ax.semilogx(k,four[i,:]/N, label=str(tplot[i]))

lg = plt.legend(loc=1)
ax.text(1.1,.475,'N = '+ str(N),fontsize=16)


ax.set_xlabel(r'Wavenumber')
ax.set_ylabel(r'$|\hat{u}|$')

plt.savefig('figs/burgers_spec_'+str(N))

