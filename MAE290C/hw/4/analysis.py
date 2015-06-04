# MAE 290C, HW 4
# Cesar B Rocha

# This code plots the magnitude of the reflection
#    coefficient for various artificial computational BC


import numpy as np
from numpy import pi, cos, sin, tan, exp
import matplotlib.pyplot as plt

plt.close('all')

kdx = np.linspace(0.001,pi/2,100)
expkx = lambda n: exp(1j*n*kdx)

# (1): homogeneous Dirichlet
r1 = np.ones_like(kdx)

# (2): Linear extrapolation
r2 = (2.*expkx(1.)-expkx(2.) - 1.) / (1. + 2.*expkx(-1.) + expkx(-2.))

# (3): Quadratic extrapolation
r3 = (3.*expkx(1.)-3.*expkx(2.)+expkx(3.)- 1.) / (1. + 3.*expkx(-1.) +\
        3.*expkx(-2.) + expkx(-3.))

# (4): homonegenous Neumann
r4 = exp(1j*(kdx-pi)/2.)*tan(kdx/2.)

# (5): first-oder, up-wind, convective bc
c = 1.
gamma = lambda A: 1/(1j*kdx) 
r5 = lambda A: (gamma(A)*expkx(1.)-gamma(A)-1.) / (1.+gamma(A) + gamma(A)*expkx(-1.)) 
r5_1 = r5(.1)
r5_3 = r5(3.)

# (5): second-order, up-wind, convective bc
c = 1.
gamma = lambda A: .5/(1j*kdx) 
r6 = lambda A: (4.*gamma(A)*expkx(1.)-3.*gamma(A)-gamma(A)*expkx(2.)-1.) /\
        (1.+3.*gamma(A) + 4.*gamma(A)*expkx(-1.)+gamma(A)*expkx(-2.))
r6_1 = r6(.1)
r6_3 = r6(3.)

## plotting
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
ax.plot(kdx,np.abs(r1),label='(1): Homogeneous Dirichlet')
ax.plot(kdx,np.abs(r2),label='(2): Linear extrapolation')
ax.plot(kdx,np.abs(r3),label='(3): Quadratic extrapolation')
ax.plot(kdx,np.abs(r4),label='(4): Homogeneous Neumann')
ax.plot(kdx,np.abs(r5_3),label='(5): First-order convective')
ax.plot(kdx,np.abs(r6_3),label='(6): Second-order convective')
plt.legend(loc=(.025,.4))
ax.set_ylim(0,1.2)
ax.set_xlim(0,pi/2)
ax.set_xlabel(r'$k \Delta x$')
ax.set_ylabel(r'$|r|$')
ax.set_xticks([0.,pi/6,pi/3,pi/2])
ax.set_xticklabels([r'$0$',r'$\pi/6$',r'$\pi/3$',r'$\pi/2$'])
plt.savefig('hw4.png')

