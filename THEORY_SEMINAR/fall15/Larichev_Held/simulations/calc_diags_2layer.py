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

# invert for ph
def qh2ph(a,qh):
    """ Inverts qh to ph using model's inversion
            object a (m.a) """
    return np.einsum("ij...,j...->i...",a,qh)

def psi_tau_h(ph):
    """ Calculates barotropic and
            baroclinic streamfunction """
    return (ph[0]+ph[1])/2., (ph[0]-ph[1])/2.

def calc_tfh(pph,pth,U,k):
    """ calculates thickness flux,
        or productions term """
    return  U*(m.rd**-2)*(1j*k*pph.conj()*pth).real

def jacob_h(model,Ah,Bh):
    """ calculates the Jacobian """

    Ax = np.fft.irfft2(1j*model.k*Ah)
    Bx = np.fft.irfft2(1j*model.k*Bh)
    Ay = np.fft.irfft2(1j*model.l*Ah)
    By = np.fft.irfft2(1j*model.l*Bh)
    J = Ax*By - Ay*Bx
    return model.filtr*np.fft.rfft2(J)

# kinetic energy
kefile = np.load(simulation+'/kinetic_energy.npz')

files =  glob.glob(pathi+"*.h5")
t, ke = [], []
Ebt, Ebc, Epe = [], [], []
E_bt_bc = []
E_ek_bt = []
k = 0

for file in files[:]:
    snap = h5py.File(file,'r')

    qh = snap['qh'][:]
    ph = qh2ph(m.a,qh)

    pph, pth = psi_tau_h(ph)

    # calculates the thickness flux
    prod = calc_tfh(pph,pth,m.U[0],m.k)

    ki, prod = diagnostic_tools.calc_ispec(m,prod)

    try:
        prod_k = np.vstack([prod_k,prod])
    except:
        prod_k = prod

    Ebt = np.hstack([Ebt, diagnostic_tools.spec_var(m,m.wv*pph)/2.])
    Ebc = np.hstack([Ebc, diagnostic_tools.spec_var(m,m.wv*pth)/2.])
    Epe = np.hstack([Epe, diagnostic_tools.spec_var(m,(1./m.rd)*pth)/2.])

    nlh = pth.conj()*jacob_h(m,pph,-m.wv2*pth)
    E_bt_bc = np.hstack([E_bt_bc, nlh.real.sum()])
    E_ek_bt =  np.hstack([E_ek_bt, diagnostic_tools.spec_var(m,m.wv*pth)])


ki, ti = np.meshgrid(ki,kefile['t'])

prod_k = np.ma.masked_array(prod_k, prod_k<=.0)

# the fractional energy production at deformation radius scales
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
fs = 20.
lw=2.
aph = .5

dir_figs = 'figs/'


plt.figure(figsize=(13,9))
plt.pcolormesh(ki,ti,prod_k/m.M**2, norm = LogNorm())
plt.xscale('log')
plt.xlim(1.,1e2)
plt.ylim(1,260)
plt.clim(1.e-7,1.e-3)
plt.xlabel('Wavenumber',fontsize=20.)
plt.xticks([1.,10.,50.,100.])
plt.ylabel('time',fontsize=20.)
plt.colorbar()
plt.savefig(dir_figs+"prod_spinup")

# energy
plt.figure(figsize=(16,9))
plt.plot(Ebt,linewidth=lw,alpha=aph,label="Barotropic")
plt.plot(Ebc+Epe,linewidth=lw,alpha=aph,label="Baroclinic")
plt.xlabel("time",fontsize=fs)
plt.ylabel("energy",fontsize=fs)
plt.legend(fontsize=fs)
plt.savefig(dir_figs+"energy_btbc")


# baroclinic energy energy
plt.figure(figsize=(16,9))
plt.plot(Ebc,linewidth=lw,alpha=aph,label="Kinetic energy")
plt.plot(Epe,linewidth=lw,alpha=aph,label="Potential energy")
plt.xlabel("time",fontsize=fs)
plt.ylabel("energy",fontsize=fs)
plt.legend(fontsize=fs)
plt.savefig(dir_figs+"energy_bc")
