import numpy as np
import matplotlib.pyplot as plt
import os
import h5py

import pyqg
import LH95_params

def init_condition(model,sig=1.e-7):
    """ White noise spectrum with amplitude sig """
    return sig*np.random.randn(model.nz,model.nx,model.ny)

def _init_save_snapshots(self,path):

    self.fno = path

    if not os.path.isdir(self.fno):
        os.makedirs(self.fno)
        os.makedirs(self.fno+"/snapshots/")

def _file_exist(self, fno):
    if os.path.exists(fno):
        if self.overwrite:
            os.remove(fno)
        else: raise IOError("File exists: {0}".format(fno))


def _save_snapshots(self, fields=['qh','u','v']):

    """ Save snapshots of fields """


    fno = self.fno + '/snapshots/{:015.0f}'.format(self.t)+'.h5'

    _file_exist(self,fno)

    h5file = h5py.File(fno, 'w')

    for field in fields:
        h5file.create_dataset(field, data=eval("self."+field))

    h5file.close()

def _save_diagnostics(self):

    """ Save diagnostics """

    fno = self.fno + 'diagnostics.h5'

    _file_exist(self,fno)

    h5file = h5py.File(fno, 'w')

    for key in self.diagnostics.keys():
        h5file.create_dataset(key, data=self.get_diagnostic(key))
      
    h5file.close()

def _save_setup(self,):

    """Save setup  """

    fno = self.fno + 'setup.h5'

    _file_exist(self,fno)

    h5file = h5py.File(fno, 'w')
    
    h5file.create_dataset("grid/nx", data=(self.nx),dtype=int)
    h5file.create_dataset("grid/ny", data=(self.ny),dtype=int)
    h5file.create_dataset("grid/nz", data=(self.nz),dtype=int)
    h5file.create_dataset("grid/x", data=(self.x))
    h5file.create_dataset("grid/y", data=(self.y))
    h5file.create_dataset("grid/H", data=(self.Hi))
    h5file.create_dataset("grid/wv", data=self.wv)

    h5file.create_dataset("basestate/U", data=self.Ubg)
    h5file.create_dataset("basestate/V", data=self.Vbg)

    h5file.create_dataset("constants/beta", data=self.f)
    h5file.create_dataset("constants/f", data=self.beta)

    h5file.close()


#
# Simulation
#

# the basic parameters
params = LH95_params.LHTwoLayer()


# model class
m = pyqg.LayeredModel(nx=params.Nx, nz=params.Nz, U = np.array([params.U1,params.U2]),
                              V=np.array([0.,0.]),L=params.L,f=params.f0,beta=params.beta,
                              rd=params.Ld, H = np.array([params.H1,params.H2]),rek=params.rek,
                              dt=0.001,tmax=params.tmax,twrite=params.twrite, 
                              tavestart=params.tavestart, ntd=4, delta=params.delta,
                              hb=None,logfile=None)

m.set_q(init_condition(m, sig=1.e-1))

m.overwrite = True

# just for testing
m.tmax = m.tmax

ke = []
t = []
patho = "REFERENCE"
k = 1

kmean = 0
_init_save_snapshots(m,patho)

for i in m.run_with_snapshots(tsnapstart=0, tsnapint=1.):

    ke.append(m._calc_ke())
    t.append(m.t)
    
    if m._calc_cfl()> .3:
        m.dt = m.dt/5.
       
    _save_snapshots(m)

_save_diagnostics(m)
_save_setup(m)

