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
#simulation = "SIZE/1_4"

#pathi = simulation+"/LONG/snapshots/"
pathi = simulation+"/snapshots/"


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

def calc_tfh(pph,pth,U,k):
    """ calculates thickness flux,
        or productions term """
    return  U*(m.rd**-2)*(1j*k*pph.conj()*pth).real

def jacob_h(model,Ah,Bh):
    """ calculates the Jacobian Jhat(A,B) """

    Ax = np.fft.irfft2(model.ik*Ah)
    Bx = np.fft.irfft2(model.ik*Bh)
    Ay = np.fft.irfft2(model.il*Ah)
    By = np.fft.irfft2(model.il*Bh)

    J = Ax*By - Ay*Bx
    return np.fft.rfft2(J)

def jacob_h2(k,l,Ah,Bh):
    """ calculates the Jacobian Jhat(A,B) """

    Ax = np.fft.ifft2(1j*k*Ah)
    Bx = np.fft.ifft2(1j*k*Bh)
    Ay = np.fft.ifft2(1j*l*Ah)
    By = np.fft.ifft2(1j*l*Bh)

    J = Ax*By - Ay*Bx

    return np.fft.fft2(J)

def ph2uv(model,ph):
    """ calculates u and v from ph """
    u = np.fft.irfft2(-model.il*ph)
    v = np.fft.irfft2( model.ik*ph)
    return u, v

def advect(u,v,q):
    uq = u*q
    vq = v*q
    return m.ik*(np.fft.rfft2(uq)) + m.il*np.fft.rfft2(vq)

# kinetic energy
#kefile = np.load(simulation+'/kinetic_energy.npz')

files =  glob.glob(pathi+"*.h5")
t, ke = [], []
Ebt, Ebc = [], []
Epe = []

Ebt_I = []
Ebt_II = []
Ebt_III = []
Ebt_IV = []
Ebt_V = []

Ebc_I = []
Ebc_IIa = []
Ebc_IIb = []
Ebc_IIIa = []
Ebc_IIIb = []
Ebc_IV = []
Ebc_V = []

k = 0

#m.il = np.flipud(m.il)

kx = np.hstack([np.arange(-m.Nx/2,0),np.arange(0,m.Nx/2)])
ky = kx.copy()

kx,ky = np.meshgrid(kx,ky)

#        self.add_diagnostic('KEflux',
#            description='spectral flux of kinetic energy',
#            function= (lambda self:
#              np.real(self.del1*self.ph[0]*np.conj(self.Jpxi[0])) +
#              np.real(self.del2*self.ph[1]*np.conj(self.Jpxi[1])) )
#        )

for file in files[30:]:

    snap = h5py.File(file,'r')

    qh = snap['qh'][:]
    ph = qh2ph(m, qh)

    pph, pth = btbc_h(ph)
    qph, qth = btbc_h(qh)

    snap.close()

    # calculate the spectrum
    try:
        Ebt +=  np.abs(m.wv*pph/(m.M))**2
        Ebc +=  np.abs(m.wv*pth/m.M)**2
        Epe +=  np.abs((1./m.rd)*pth/m.M)**2
    except:
        Ebt =  np.abs(m.wv*pph/m.M)**2
        Ebc =  np.abs(m.wv*pth/m.M)**2
        Epe =  np.abs((1./m.rd)*pth/m.M)**2

    # the RHS of the barotropic spectrum equation
    #Jpph = jacob_h(m,pph,-m.wv2*pph)
    pp, pt = np.fft.irfft2(pph), np.fft.irfft2(pth)
    qp, qt = np.fft.irfft2(qph), np.fft.irfft2(qth)

    #pph2, pth2 = np.fft.fft2(pp), np.fft.fft2(pt)
    #qph2, qth2 = np.fft.fft2(qp), np.fft.fft2(qt)

    up, vp = ph2uv(m, pph)
    qp = np.fft.irfft2(-m.wv2*pph)
    Jpph = advect(up,vp,qp)
    Jpphc = Jpph.conj()
    pphc = pph.conj()

    #Jtth = jacob_h(m,pth,-m.wv2*pth)
    ut, vt = ph2uv(m, pth)
    qt = np.fft.irfft2(-m.wv2*pth)
    Jtth = advect(ut,vt,qt)
    Jtthc = Jtth.conj()
    pthc = pth.conj()

    #Jpph2 = jacob_h2(kx,ky,pph2,qph2)
    #Jpph2 = Jpph2.conj()

    Jpth = jacob_h(m,pph,-m.wv2*pth)
    Jpthc = Jpth.conj()

    Jtph = jacob_h(m,pth,qph)
    Jtphc = Jtph.conj()

    try:

        Ebt_I   += (pph*Jpphc).real/m.M**2
        #Ebt_I2  += (pph2.conj()*Jpph2).real/m.M**2

        Ebt_II  += (pph*Jtthc).real/m.M**2
        Ebt_III += -m.U[0]*m.k*m.wv2*( (1j*pphc*pth).real )/m.M**2
        Ebt_V   += -.5*m.rek*m.wv2*(pphc*ph[-1] ).real/m.M**2

        Ebc_I +=    (pth*Jtphc).real/m.M**2
        Ebc_IIa +=  (pth*Jpthc).real/m.M**2
        Ebc_IIb +=   -(m.rd**-2)*(pthc*jacob_h(m,pph,pth)).real/m.M**2
        Ebc_IIIa +=  -m.U[0]*m.wv2*m.k*(1j*pthc*pph).real/m.M**2
        Ebc_IIIb +=  m.U[0]*(m.rd**-2)*m.k*(1j*pthc*pph).real/m.M**2

        Ebc_V    +=  .5*m.rek*m.wv2*( (pthc*ph[-1]).real)/m.M**2

    except:

        Ebt_I   = (pph*Jpphc).real/m.M**2
        #Ebt_I2  = (pph2.conj()*Jpph2).real/m.M**2

        Ebt_II  = (pph*Jtthc).real/m.M**2
        Ebt_III = -m.U[0]*m.k*m.wv2*( (1j*pphc*pth).real )/m.M**2
        Ebt_V   = -.5*m.rek*m.wv2*(pphc*ph[-1] ).real/m.M**2

        Ebc_I   =  (pth*Jtphc).real/m.M**2
        Ebc_IIa =  (pth*Jpthc).real/m.M**2
        Ebc_IIb =   -(m.rd**-2)*(pthc*jacob_h(m,pph,pth)).real/m.M**2
        Ebc_IIIa =  -m.U[0]*m.wv2*m.k*(1j*pth*pphc).real/m.M**2
        Ebc_IIIb =  m.U[0]*(m.rd**-2)*m.k*(1j*pthc*pph).real/m.M**2

        Ebc_V    =  .5*m.rek*m.wv2*( (pthc*ph[-1]).real)/m.M**2

    k+=1

Ebt = Ebt/k
Ebc = Ebc/k
Epe = Epe/k
Ebt_I = Ebt_I/k
#Ebt_Ib = Ebt_Ib/k
Ebt_II = Ebt_II/k
Ebt_III = Ebt_III/k
Ebt_V = Ebt_V/k

Ebc_I = Ebc_I/k
Ebc_IIa = Ebc_IIa/k
Ebc_IIb = Ebc_IIb/k
Ebc_IIIa = Ebc_IIIa/k
Ebc_IIIb = Ebc_IIIb/k

Ebc_V = Ebc_V/k

Ebt =   Ebt.sum(axis=0)
Ebc =   Ebc.sum(axis=0)
Epe =   Epe.sum(axis=0)
# Ebt_I = Ebt_I.sum(axis=0)
# Ebt_II = Ebt_II.sum(axis=0)
# Ebt_III = Ebt_III.sum(axis=0)
# Ebt_V = Ebt_V.sum(axis=0)

#ki, Ebt =  diagnostic_tools.calc_ispec(m, Ebt)
#_, Ebc =   diagnostic_tools.calc_ispec(m, Ebc)
#_, Epe =   diagnostic_tools.calc_ispec(m, Epe)
#ki, Ebt_I = diagnostic_tools.calc_ispec(m, Ebt_I)
#_, Ebt_Ib = diagnostic_tools.calc_ispec(m, Ebt_Ib)
#_, Ebt_II = diagnostic_tools.calc_ispec(m, Ebt_II)
#_, Ebt_III = diagnostic_tools.calc_ispec(m, Ebt_III)
#_, Ebt_V = diagnostic_tools.calc_ispec(m, Ebt_V)

# _, Ebc_I = diagnostic_tools.calc_ispec(m, Ebc_I)
# _, Ebc_IIa = diagnostic_tools.calc_ispec(m, Ebc_IIa)
# _, Ebc_IIb = diagnostic_tools.calc_ispec(m, Ebc_IIb)
# _, Ebc_IIIa = diagnostic_tools.calc_ispec(m, Ebc_IIIa)
# _, Ebc_IIIb = diagnostic_tools.calc_ispec(m, Ebc_IIIb)
# _, Ebc_V = diagnostic_tools.calc_ispec(m, Ebc_V)

ki = m.kk

Ebt_I = Ebt_I.sum(axis=0)
Ebt_II = Ebt_II.sum(axis=0)
Ebt_III = Ebt_III.sum(axis=0)
Ebt_V = Ebt_V.sum(axis=0)

Ebc_I = Ebc_I.sum(axis=0)
Ebc_IIa = Ebc_IIa.sum(axis=0)
Ebc_IIb = Ebc_IIb.sum(axis=0)
Ebc_IIIa = Ebc_IIIa.sum(axis=0)
Ebc_IIIb = Ebc_IIIb.sum(axis=0)
Ebc_V = Ebc_V.sum(axis=0)

Ebt_all = [Ebt_I, Ebt_II, Ebt_III, Ebt_V ]
#Ebt_all.append(np.vstack(Ebt_all).sum(axis=0))

Ebc_all = [ Ebc_I, Ebc_IIa, Ebc_IIb, Ebc_IIIa, Ebc_IIIb, Ebt_V]
#Ebc_all.append(np.vstack(Ebc_all).sum(axis=0))

plt.figure(figsize=(12,8))
labels = ['I','II','III','V']
[plt.semilogx(ki, term,linewidth=2.,alpha=.5) for term in  Ebt_all]
plt.legend(labels, loc='upper right',fontsize=20)
plt.xlim([m.kk.min(), m.kk.max()])
plt.xlabel(r'k ',fontsize=20); plt.grid()
plt.ylabel(r'Spectral density',fontsize=20)
plt.title('Barotropic Spectral Energy Transfers');
plt.savefig('barotropic_transfer')


plt.figure(figsize=(12,8))
labels = ['I','IIa','IIb','IIIa','IIIb','V']
[plt.semilogx(ki, term, linewidth=2.,alpha=.5) for term in  Ebc_all]
plt.legend(labels, loc='upper right',fontsize=20)
plt.xlim([m.kk.min(), m.kk.max()])
plt.xlabel(r'k ',fontsize=20); plt.grid()
plt.ylabel(r'Spectral  density',fontsize=20)
plt.title('Baroclinic Spectral Energy Transfers',fontsize=20);
plt.savefig('baroclinic_transfer')

# energy flux
ep_bt = np.cumsum(Ebt_I)
ep_bc = np.cumsum(Ebc_IIb+Ebc_IIIa)


plt.figure(figsize=(12,8))
plt.semilogx(ki,ep_bt,linewidth=2.,alpha=.5,label="Batropic")
plt.semilogx(ki,ep_bc,linewidth=2.,alpha=.5,label="Baroclinic")
plt.xlabel(r'k ',fontsize=20); plt.grid()
plt.ylabel(r'Spectral  energy flux',fontsize=20)
plt.ylim(-.02,.02)
plt.xlim(1.,m.kk.max())
plt.savefig('energy_flux')


# the residual
plt.figure()
plt.semilogx(ki,Ebt_all[-1]+Ebc_all[-1],label="Resid.")
plt.legend(loc='upper right')
plt.xlim([m.kk.min(), m.kk.max()])
plt.title('Spectral Energy Transfers');
plt.ylabel(r'Spectral densities',fontsize=20)

plt.savefig('spec_bc')

#Ebt = Ebt.mean(axis=0)
#Ebc = Ebc.mean(axis=0)
#Ebc2 = Ebc2.mean(axis=0)
#Epe = Epe.mean(axis=0)

# plot spectrum

# plot
aph = .5
ks = np.array([1.,50.])

plt.figure(figsize=(12,8))
plt.loglog(m.kk,Ebt,label=r'$E_\psi(k)$',linewidth=3.,alpha=aph)
plt.loglog(m.kk,Ebc+Epe,label=r'$E_\tau(k)$',linewidth=3.,alpha=aph)
plt.loglog(ks,ks**(-5/3.)/3.2,'k--',linewidth=2.)
plt.text(1.3,.3,r'$k^{-5/3}$',fontsize=20)
plt.xlim([m.kk.min(), m.kk.max()])
plt.xlabel(r'$k$',fontsize=20); plt.grid()
plt.ylabel('Spectral energy density',fontsize=20);
plt.legend(fontsize=20)
plt.ylim(1.e-6,1.)
plt.savefig('spec')
