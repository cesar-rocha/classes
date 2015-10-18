import numpy as np
from numpy import pi
from pyearth import earth

# parameters
class LHTwoLayer():
    """ Standard parameters for ACC-like
            2-layer simulation """

    def __init__(self):

        self.Nx = 256          # number of grid points
        self.Nz = 2            # number of layers
        self.L =  2*pi         # length scale of box    [m]
        self.Ld = 1./50.       # deformation scale      [m]
        self.kd = 1./self.Ld   # deformation wavenumber [m^-1]
        self.H1 = 0.5          # upper layer thickness  [m]
        self.H2 = 0.5          # lower layer thickness  [m]
        self.delta = self.H1/self.H2

        self.b = 10*self.Ld    # deccay scale of topographic feature [m]
        self.A = .0            # amplitude of topographic feature h/H2 [unitless]

        self.U1 = .005          # upper layer zonal velocity [m/s]
        self.U2 = -.005         # lower layer zonal velicity [m/s]
        self.Us = self.U1-self.U2  # vertical shear             [m/s]

        self.rek = .04          # linear bottom drag coeff.  [s^-1]

        self.lat = 45.          # latitude

        betaplane = earth.Betaplane(lat=self.lat)
        self.f0  = betaplane.f0
        self.beta = 0.
         
        self.Ti = self.Ld/(abs(self.Us))  # estimate of most unstable e-folding time scale [s]

        self.dx = self.L/(self.Nx-1)            # resolution
        self.x = np.linspace(0.,self.L,self.Nx)
        self.y = self.x.copy()
        self.x,self.y = np.meshgrid(self.x,self.y)

        # topography
        self.eta = self.A*np.exp( -((self.x-self.x[self.Nx/2,self.Nx/2])/self.b)**2 
                     -((self.y-self.y[self.Nx/2,self.Nx/2])/self.b)**2 )

        self.dt = self.Ti/1000.      # time-step [s]
        self.tmax = 4000             # simulation time [s]
        self.twrite=200
        self.tavestart=1500.

