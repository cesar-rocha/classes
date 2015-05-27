import numpy as np
import matplotlib.pyplot as plt
import spectral_model as model

reload(model)

# the model object
m = model.BTModel(Lx=2.*np.pi,nx=256, tmax = 50,dt = 0.005, \
        use_fftw=True, ntd=4,use_filter=True,nu=0.)


#m.set_q(np.cos(0.*m.x + 10.*m.y))
# McWilliams 84 IC condition
#fk = m.kappa != 0
#ckappa = np.zeros_like(m.kappa2)
#ckappa[fk] = np.sqrt( m.kappa2[fk]*(1. + (m.kappa2[fk]/36.)**2) )**-1
#
#nhx,nhy = m.kappa2.shape
#
#Pi_hat = np.random.randn(nhx,nhy)*ckappa +1j*np.random.randn(nhx,nhy)*ckappa
#
#Pi = m.ifft2( Pi_hat[:,:] )
#Pi = Pi - Pi.mean()
#Pi_hat = m.fft2( Pi )
#KEaux = m.spec_var(m.kappa*Pi_hat )
#
#pih = ( Pi_hat/np.sqrt(KEaux) )
#qih = -m.kappa2*pih
#qi = m.ifft2(qih)
#m.set_q(qi)


m.run()

#m.run_with_snapshots(tsnapstart=0., tsnapint=50)

# run the model and plot some figs
#plt.rcParams['image.cmap'] = 'RdBu'
#
#plt.ion()
#
#for snapshot in m.run_with_snapshots(tsnapstart=0, tsnapint=100*m.dt):
#
#    plt.clf()
#    p1 = plt.contourf(m.q,np.linspace(-30,30,12))
#    plt.clim([-30., 30.])
#    plt.title('t='+str(m.t))
#
#    plt.xticks([])
#    plt.yticks([])
#
#    plt.pause(0.01)
#
#    plt.draw()
#
#plt.show()
#plt.ion()
