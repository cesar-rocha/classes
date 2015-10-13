import sys
import glob
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import pyqg
from pyqg import diagnostic_tools
import LH95_params

plt.close('all')

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

# invert for ph
def qh2ph(a,qh):
    """ Inverts qh to ph using model's inversion
            object a (m.a) """
    return np.einsum("ij...,j...->i...",a,qh)

def btbc_h(ph):
    """ Calculates barotropic and 
            baroclinic streamfunction """
    return (ph[0]+ph[1])/2., (ph[0]-ph[1])/2.

def calc_tfh(pph,pth,U,k):
    """ calculates thickness flux,
        or productions term """
    return  U*(m.rd**-2)*(1j*k*pph.conj()*pth).real

def jacob_h(model,Ah,Bh):
    """ calculates the Jacobian """

    Ax = np.fft.irfft2(model.filtr*(1j*model.k*Ah))
    Bx = np.fft.irfft2(model.filtr*(1j*model.k*Bh))
    Ay = np.fft.irfft2(model.filtr*(1j*model.l*Ah))
    By = np.fft.irfft2(model.filtr*(1j*model.l*Bh))
    J = Ax*By - Ay*Bx
    return model.filtr*np.fft.rfft2(J)

# kinetic energy
kefile = np.load(simulation+'/kinetic_energy.npz')

files =  glob.glob(pathi+"*.h5")
t, ke = [], []
Ebt, Ebc = [], []
Epe = []

Ebt_I = []
Ebt_II = []
Ebt_III = []
Ebt_IV = []
Ebt_V = []

k = 0

for file in files[500:]:
    snap = h5py.File(file,'r')

    qh = snap['qh'][:]
    ph = qh2ph(m.a,qh)

    pph, pth = btbc_h(ph)
    qph, qth = btbc_h(qh)

    # calculate the spectrum
    try:
        Ebt +=  np.abs(m.wv*pph/(m.M))**2
        Ebc +=  np.abs(m.wv*pth/m.M)**2
        Epe  += np.abs((1./m.rd)*pth/m.M)**2
    except:
        Ebt =  np.abs(m.wv*pph/m.M)**2
        Ebc =  np.abs(m.wv*pth/m.M)**2
        Epe = np.abs((1./m.rd)*pth/m.M)**2

    # the RHS of the barotropic spectrum equation
    try:
        Ebt_I +=  (pph.conj()*jacob_h(m,pph,qph)).real
        Ebt_II += (pph.conj()*jacob_h(m,pth,qth)).real
        Ebt_III += -m.U[0]*m.k*m.wv2*( (1j*pph.conj()*pth).real )
        Ebt_V += -m.rek*m.wv2*( (pph.conj()*ph[-1] ).real )
    except:
        Ebt_I =   (pph.conj()*jacob_h(m,pph,qph)).real
        Ebt_II =  (pph.conj()*jacob_h(m,pth,qth)).real
        Ebt_III = -m.U[0]*m.k*m.wv2*( (1j*pph.conj()*pth).real )
        Ebt_V = -m.rek*m.wv2*( (pph.conj()*ph[-1] ).real )

    k+=1

Ebt = Ebt/k
Ebc = Ebc/k
Epe = Epe/k
Ebt_I = Ebt_I/k
Ebt_II = Ebt_II/k
Ebt_III = Ebt_III/k
Ebt_V = Ebt_V/k

ki, Ebt =  diagnostic_tools.calc_ispec(m, Ebt) 
_, Ebc =   diagnostic_tools.calc_ispec(m, Ebc) 
_, Epe =   diagnostic_tools.calc_ispec(m, Epe) 
_, Ebt_I = diagnostic_tools.calc_ispec(m, Ebt_I) 
_, Ebt_II = diagnostic_tools.calc_ispec(m, Ebt_II) 
_, Ebt_III = diagnostic_tools.calc_ispec(m, Ebt_III) 
_, Ebt_V = diagnostic_tools.calc_ispec(m, Ebt_V) 


#Ebt = Ebt.mean(axis=0)
#Ebc = Ebc.mean(axis=0)
#Ebc2 = Ebc2.mean(axis=0)
#Epe = Epe.mean(axis=0)
