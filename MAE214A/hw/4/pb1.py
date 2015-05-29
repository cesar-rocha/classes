import numpy as np
import matplotlib.pyplot as plt
import scipy.io as io

plt.close('all')

# reading data
pb1 = io.loadmat('mae214A_SPR15_hw4/pb1.mat')

# (a) and (b) u at grid point (65,iz)
def umeans(iz):
    uy = pb1['u'][:,:,iz].mean(axis=1)
    ut = pb1['u'][:,65,iz]
    uyt = pb1['u'][:,:,iz].mean().repeat(uy.size)
    return ut,uy,uyt

# plotting
for iz in [2,8,20]:

    ut,uy,uyt = umeans(iz)

    fig = plt.figure(figsize=(12,6))
    plt.plot(ut,label=r'$u(t)$')
    plt.plot(uy,label=r'$<u>_y(t)$')
    plt.plot(uyt,label=r'$<u>_{yt}$')
    lg = plt.legend(loc=2)
    plt.xlabel('Time')
    plt.ylabel('u')
    plt.savefig('pb1ab_'+str(iz))

# (c) <u+>_yt
up = pb1['u'].mean(axis=(0,1))[:64]
zp = pb1['zplus']

fig = plt.figure(figsize=(8,12))
plt.semilogy(up,zp)
#plt.semilogy(zp,zp,'k--')
#plt.xlim(0.,1e3)
#plt.ylim(0,20.)
plt.xlabel(r'$z^+$')
plt.ylabel(r'$u^+$')
plt.savefig('pb1c')

# (d)
up_data = np.loadtxt('mae214A_SPR15_hw4/Uplus.dat')
fig = plt.figure(figsize=(8,12))
plt.semilogy(up_data[:,1],up_data[:,0])
#plt.loglog(zp,zp,'k--')
