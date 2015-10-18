import sys
import glob
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import pyqg
from pyqg import diagnostic_tools
import LH95_params

# the basic parameters
params = LH95_params.LHTwoLayer()

# model class
m = pyqg.LayeredModel(nx=params.Nx, nz=params.Nz, U = np.array([params.U1,params.U2]),
                              V=np.array([0.,0.]),L=params.L,f=params.f0,beta=params.beta,
                              rd=params.Ld, H = np.array([params.H1,params.H2]),rek=params.rek,
                              dt=0.001,tmax=params.tmax,twrite=params.twrite,
                              tavestart=params.tavestart, ntd=2, delta=params.delta,
                              hb=params.eta*params.H2,logfile=None)

plt.close('all')

# the basic parameters
params = LH95_params.LHTwoLayer()

# invert for ph
def qh2ph(model,qh):
    """ Inverts qh to ph using model's inversion
            object a (m.a) """
    ph = np.einsum("ij...,j...->i...",model.a,qh)

    return  ph

def btbc_h(p):
    """ Calculates barotropic and
            baroclinic streamfunction """
    return (p[0]+p[1])/2., (p[0]-p[1])/2.

def fft2(A):
    return np.fft.rfft2(A)

def ifft2(Ah):
    return np.fft.irfft2(Ah)

# the basic parameters
params = LH95_params.LHTwoLayer()

# model class
m = pyqg.LayeredModel(nx=params.Nx, nz=params.Nz, U = np.array([params.U1,params.U2]),
                              V=np.array([0.,0.]),L=params.L,f=params.f0,beta=params.beta,
                              rd=params.Ld, H = np.array([params.H1,params.H2]),rek=params.rek,
                              dt=0.001,tmax=params.tmax,twrite=params.twrite,
                              tavestart=params.tavestart, ntd=4, delta=params.delta,
                              hb=None,logfile=None)


z = np.array([1.5, .5])

# some data
pathi = "REFERENCE/snapshots/"

#snap = h5py.File(pathi+"000000000001728.h5","r")

simulation = "REFERENCE"
#simulation = "SIZE/1_4"

files =  glob.glob(pathi+"*.h5")

file = files[3000]
snap = h5py.File(file,'r')

qh = snap['qh'][:]
ph = qh2ph(m, qh)

pph, pth = btbc_h(ph)
qph, qth = btbc_h(qh)

# now transform to physical space
q = ifft2(qh)
q += m.Qy[:,np.newaxis,np.newaxis]*m.y

qp, qt = ifft2(qph), ifft2(qth) - m.U[0]*m.y
pp, pt = ifft2(pph), ifft2(pth)

plt.figure()
plt.contourf(m.x/m.rd,m.y/m.rd,qp,100,cmap="RdBu_r")
plt.xlabel(r'$x\, k_d$',fontsize=20)
plt.ylabel(r'$y\, k_d$',fontsize=20)
plt.clim([-50,50.])
plt.xlim(0,314)
plt.ylim(0,314)
plt.savefig('qbt')

plt.figure()
plt.contourf(m.x/m.rd,m.y/m.rd,qt,100,cmap="RdBu_r")
plt.xlabel(r'$x\, k_d$',fontsize=20)
plt.ylabel(r'$y\, k_d$',fontsize=20)
plt.clim([-90,90.])
plt.xlim(0,314)
plt.ylim(0,314)
plt.savefig('qbc')



plt.figure()
plt.contourf(m.x/m.rd,m.y/m.rd,pp,100,cmap="RdBu_r")
plt.xlabel(r'$x\, k_d$',fontsize=20)
plt.ylabel(r'$y\, k_d$',fontsize=20)
plt.clim([-2.,2])
plt.xlim(0,314)
plt.ylim(0,314)
plt.savefig('pbt')

plt.figure()
plt.contourf(m.x/m.rd,m.y/m.rd,pt,100,cmap="RdBu_r")
plt.xlabel(r'$x\, k_d$',fontsize=20)
plt.ylabel(r'$y\, k_d$',fontsize=20)
plt.clim([-.035,.035])
plt.xlim(0,314)
plt.ylim(0,314)
plt.savefig('pbc')
