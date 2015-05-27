from __future__ import division
import numpy as np
from numpy import pi, exp, sqrt, cos, sin

try:   
    import mkl
    np.use_fastnumpy = True
except ImportError:
    pass

try:
    import pyfftw
    pyfftw.interfaces.cache.enable() 
except ImportError:
    pass

class BTModel(object):

    """ A class that represents the 2D model """

    def __init__(

            self,
            nx = 256,
            ny=None,
            Lx=2*pi,                    # domain size
            Ly=None,                    # width
            # physical parameters
            nu = 0.,
            # timestepping parameters
            dt=.0025,               # numerical timestep
            twrite=100,             # interval for cfl and ke printout (in timesteps)
            tmax=100.,              # total time of integration
            filt=True,              # spectral filter flag
            use_fftw=True,
            ntd = 1,
            use_filter=True):

        if ny is None: ny = nx
        if Ly is None: Ly = Lx
       
        # initilize parameters

        # domain
        self.nx = nx
        self.ny = ny
        self.Lx = Lx
        self.Ly = Ly
        self.dx = Lx/nx
        self.dy = Ly/ny
        self.x = np.arange(0.,Lx,self.dx)
        self.y = np.arange(0.,Ly,self.dy)
        self.x,self.y = np.meshgrid(self.x,self.y)

        # physical
        self.nu = nu

        # time
        self.dt = dt
        self.twrite = twrite
        self.tmax = tmax
        self.t = 0.
        self.ndt = 0

        # fourier settings
        self._init_kxky()
        self.kappa2 = self.k**2 + self.l**2
        self.kappa = sqrt(self.kappa2)
        self.fnz = self.kappa2 != 0
        self.kappa2i = np.zeros_like(self.kappa2)   # inversion not defined at kappa=0
        self.kappa2i[self.fnz] = self.kappa2[self.fnz]**-1
        
        # a filter
        if use_filter:
            self.filt = spec_filt(self,cphi=0.65*pi)
        else:
            # if not use exponential filter,
            #   then dealias using 2/3 rule
            self.filt = np.ones_like(self.kappa2)
            self.filt[self.nx/3:2*self.nx/3,:] = 0.
            self.filt[:,self.ny/3:] = 0.

        # other
        self.use_fftw = use_fftw
        self.ntd = ntd

        # allocate variables
        self._allocate_variables()

        # DFT
        self._initialize_fft()

        # initialize step forward
        self._init_rk3w()

        # initialize vorticity
        self.set_q(np.random.randn(self.ny,self.nx))

    def run(self):
        """ step forward until tmax """

        while(self.t < self.tmax): 
            
            self._stepforward()

            if (self.ndt%self.twrite == 0.):
                self._printout()
        
            self.t += self.dt
            self.ndt += 1

    def run_with_snapshots(self, tsnapstart=0., tsnapint=432000.):
        """ Run the model forward until the next snapshot, then yield."""
        
        tsnapints = np.ceil(tsnapint/self.dt)
        nt = np.ceil(np.floor((self.tmax-tsnapstart)/self.dt+1)/tsnapints)
        
        while(self.t < self.tmax):
            self._stepforward()
            if (self.ndt%self.twrite == 0.):
                self._printout()
            if self.t>=tsnapstart and (self.ndt%tsnapints)==0:
                yield self.t    
            self.t += self.dt
            self.ndt += 1
        return

    def _stepforward(self):

        """ march the system forward using a RK3W-theta scheme """
 
        self.nl1h = -self.jacobian()
        self.qh = (self.L1*self.qh + self.c1*self.dt*self.nl1h).copy()
        self.qh = self.filt*self.qh

        self.nl2h = self.nl1h.copy()
        self.nl1h = -self.jacobian()
        self.qh = (self.L2*self.qh + self.c2*self.dt*self.nl1h +\
                self.d1*self.dt*self.nl2h).copy()
        self.qh = self.filt*self.qh

        self.nl2h = self.nl1h.copy()
        self.nl1h = -self.jacobian()
        self.qh = (self.L3*self.qh + self.c3*self.dt*self.nl1h +\
                self.d2*self.dt*self.nl2h).copy()
        self.qh = self.filt*self.qh

        # AB2
        #if self.ndt == 0:
        #    self.nl1h = -self.dt*self.jacobian()
        #    self.qh = self.filt*(self.qh + self.nl1h)
        #    self.nl2h = self.nl1h.copy()
        #else:
        #    self.nl1h = -self.dt*self.jacobian()
        #    self.qh = self.filt*(self.qh + 1.5* self.nl1h  -0.5 * self.nl2h)
        #    self.nl2h = self.nl1h.copy()

        # forward euler
        #self.qh = self.qh - self.dt*self.jacobian()
 
    def _allocate_variables(self):
        """ Allocate variables in memory """

        dtype_real = np.dtype('float64')            
        dtype_cplx = np.dtype('complex128')
        shape_real = (self.ny, self.nx)
        shape_cplx = (self.ny, self.nx/2+1)
            
        # vorticity
        self.q  = np.zeros(shape_real, dtype_real)
        self.qh = np.zeros(shape_cplx, dtype_cplx)
        # streamfunction
        self.p  = np.zeros(shape_real, dtype_real)
        self.ph = np.zeros(shape_cplx, dtype_cplx)
        # velocity
        self.u = np.zeros(shape_real, dtype_real)
        self.v = np.zeros(shape_real, dtype_real)
        # nonlinear-term
        self.nl1h = np.zeros(shape_cplx, dtype_cplx)
        self.nl2h = np.zeros(shape_cplx, dtype_cplx)

    def _initialize_fft(self):
        # set up fft functions for use later
        if self.use_fftw:
            
            self.fft2 = (lambda x :
                    pyfftw.interfaces.numpy_fft.rfft2(x, threads=self.ntd,\
                            planner_effort='FFTW_ESTIMATE'))
            self.ifft2 = (lambda x :
                    pyfftw.interfaces.numpy_fft.irfft2(x, threads=self.ntd,\
                            planner_effort='FFTW_ESTIMATE'))
        else:

            self.fft2 =  (lambda x : np.fft.rfft2(x))
            self.ifft2 = (lambda x : np.fft.irfft2(x))

    def _init_kxky(self):
        """ Calculate wavenumbers """

        self.dl = 2.*pi/self.Ly
        self.dk = 2.*pi/self.Lx
        self.ll = self.dl*np.append( np.arange(0.,self.ny/2),
                np.arange(-self.ny/2,0.) )
        self.kk = self.dk*np.arange(0.,self.nx/2+1)
        self.k,self.l = np.meshgrid(self.kk,self.ll)
        self.kj = 1j*self.k
        self.lj = 1j*self.l

    def _invert(self):
        """ Compute streamfunction from vorticity """
        self.ph = -self.kappa2i*self.qh

    def set_q(self,q):
        """ Initialize vorticity """ 

        self.q = q
        self.qh = self.fft2(self.q)
        self._invert()
        self.ph = self.filt * self.ph
        self.ph = self.ph/(sqrt(2.*self._calc_ke()))
        self.qh = -self.kappa2*self.ph 

    def _printout(self):
        """Output some basic stats."""
        ke = self._calc_ke()
        cfl = self._calc_cfl()
        print 't=%16d, cfl=%5.6f, ke=%9.9f' % (
               self.t, cfl, ke)
        assert cfl<1., "CFL condition violated"

    def jacobian(self):

        """ compute the jacobian in conservative form """

        self._invert()


        # dealias

        self.q = self.ifft2(self.qh)
        self.u = self.ifft2(-self.lj*self.ph) 
        self.v = self.ifft2( self.kj*self.ph)

        jach = self.kj*self.fft2(self.u*self.q) +\
                self.lj*self.fft2(self.v*self.q)

        # dealias

        return jach

    # step forward
    def _init_rk3w(self):
        """ Initialize stuff for marching scheme """

        self.a1, self.a2, self.a3 = 29./96., -3./40., 1./6.
        self.b1, self.b2, self.b3 = 37./160., 5./24., 1./6.
        self.c1, self.c2, self.c3 = 8./15., 5./12., 3./4.
        self.d1, self.d2 = -17./60., -5./12.

        self.Lin = -self.nu*self.kappa2*self.dt 
        self.L1 = ( (1. + self.a1*self.Lin)/(1. - self.b1*self.Lin) )     
        self.L2 = ( (1. + self.a2*self.Lin)/(1. - self.b2*self.Lin) )
        self.L3 = ( (1. + self.a2*self.Lin)/(1. - self.b3*self.Lin) )

    # some diagnostics
    def _calc_cfl(self):
        return np.abs(
            np.hstack([self.u, self.v])).max()*self.dt/self.dx

    def _calc_ke(self):
        ke = .5*self.spec_var(self.kappa*self.ph)
        return ke.sum()

    def _calc_ens(self):
        ens = .5*self.spec_var(self.kappa2*self.ph)

    def _calc_eddy_time(self):
        """ estimate the eddy turn-over time in days """

        ens = .5 * self.spec_var(self.kappa2*self.ph)

        return 2.*pi*np.sqrt( ens**-1 )

    def spec_var(self,ph):
        """ compute variance of p from Fourier coefficients ph """
        var_dens = 2. * np.abs(ph)**2 / (self.nx*self.ny)**2 
        # only half of coefs [0] and [nx/2+1] due to symmetry in real fft2
        var_dens[:,0],var_dens[:,-1] = var_dens[:,0]/2.,var_dens[:,-1]/2.
        return var_dens.sum()


#### OLD STUFF ###

# some off-class diagnostics 
#def spec_var(self,ph):
#    """ compute variance of p from Fourier coefficients ph """
#    var_dens = 2. * np.abs(ph)**2 / (self.nx*self.ny)**2 
#    # only half of coefs [0] and [nx/2+1] due to symmetry in real fft2
#    var_dens[:,0],var_dens[:,-1] = var_dens[:,0]/2.,var_dens[:,-1]/2.
#    return var_dens.sum()

# spectral filter
def spec_filt(self,cphi=0.65*pi):
    """ Set spectral filter """
    wvx=sqrt((self.k*self.dx)**2.+(self.l*self.dy)**2.)
    filtr = exp(-23.6*(wvx-cphi)**4.)  
    filtr[wvx<=cphi] = 1.                   
    return filtr

