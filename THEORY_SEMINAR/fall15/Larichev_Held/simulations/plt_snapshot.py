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

# some data
pathi = "REFERENCE/snapshots/"

#snap = h5py.File(pathi+"000000000001728.h5","r")

simulation = "REFERENCE"
pathi = simulation+"/snapshots/"


def q_bcbt_h(qh):
    """ Calculates barotropic and 
            baroclinic  pv """
    return (qh[0]+qh[1])/2., (qh[0]-qh[1])/2.


file = 'REFERENCE/snapshots/000000000003668.h5'
snap = h5py.File(file,'r')

qh = snap['qh'][:]

qbth, qbch = q_bcbt_h(qh)

qbt, qbc = np.fft.irfft2(qbth), np.fft.irfft2(qbch)

plt.figure()
plt.pcolormesh(m.x/m.rd,m.y/m.rd,qbc)
plt.xlabel(r'$x\, m_1$')
plt.ylabel(r'$y\, m_1$')
plt.clim([-100,100.])
plt.xlim(0,310)
plt.ylim(0,310)


plt.figure()
plt.pcolormesh(m.x/m.rd,m.y/m.rd,qbt)
plt.xlabel(r'$x\, m_1$')
plt.ylabel(r'$y\, m_1$')
plt.clim([-40,40.])
plt.xlim(0,310)
plt.ylim(0,310)

