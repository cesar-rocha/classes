import numpy as np
import matplotlib.pyplot as plt
import scipy.io as io

plt.close('all')

path = 'mae214A_SPR15_hw4/'
reystress = np.loadtxt(path+'chan395.reystress')
means = np.loadtxt(path+'chan395.means')
re2000 = np.loadtxt(path+'Re2000.stat')
myres = np.load('myrestress.npz')

# reystress
# y,  y+ , R_uu, R_vv, R_ww, R_uv, R_uw , R_vw

# means
# y, y+, Umean, dUmean/dy, Wmean, dWmean/dy, Pmean



fig = plt.figure(figsize=(12,8))
plt.plot(myres['zp'],myres['uu'][:64],label=r'$<uu>$')
plt.plot(myres['zp'],myres['ww'][:64],label=r'$<vv>$')
plt.plot(myres['zp'],myres['vv'][:64],label=r'$<ww>$')
plt.plot(myres['zp'],myres['uw'][:64],label=r'$<uv>$')
plt.plot(myres['zp'],myres['uv'][:64],label=r'$<uw>$')
plt.plot(myres['zp'],myres['vw'][:64],label=r'$<vw>$')
plt.plot(myres['zp'],myres['prod'][:64],'--',label=r'$\mathcal{P}$') # check sign
#plt.xlim(0.,1e3)
plt.ylim(-2.,8.)
lg = plt.legend(loc=1)
plt.xlabel(r'$y^+$')
plt.ylabel(r'Reynolds stresses')
plt.savefig('pb2a_my')


fig = plt.figure(figsize=(12,8))
plt.plot(reystress[:,1],reystress[:,2],label=r'$<uu>$')
plt.plot(reystress[:,1],reystress[:,3],label=r'$<vv>$')
plt.plot(reystress[:,1],reystress[:,4],label=r'$<ww>$')
plt.plot(reystress[:,1],reystress[:,5],label=r'$<uv>$')
plt.plot(reystress[:,1],reystress[:,6],label=r'$<uw>$')
plt.plot(reystress[:,1],reystress[:,7],label=r'$<vw>$')
plt.plot(means[:,1],means[:,6],'--',label=r'$\mathcal{P}$') # check sign
#plt.xlim(0.,1e3)
plt.ylim(-2.,8.)
lg = plt.legend(loc=1)
plt.xlabel(r'$y^+$')
plt.ylabel(r'Reynolds stresses')
plt.savefig('pb2a')

# Re2000
# y/h, y+, U+, u'+, v'+, w'+, -Om_z+, om_x'+, om_y'+, om_z'+, uv'+, uw'+, vw'+, pr'+, ps'+, psto'+, p'
fig = plt.figure(figsize=(12,8))
plt.plot(re2000[:,1],re2000[:,3],label=r'$<uu>$')
plt.plot(re2000[:,1],re2000[:,4],label=r'$<vv>$')
plt.plot(re2000[:,1],re2000[:,5],label=r'$<ww>$')
plt.plot(re2000[:,1],re2000[:,10],label=r'$<uv>$')
plt.plot(re2000[:,1],re2000[:,11],label=r'$<uw>$')
plt.plot(re2000[:,1],re2000[:,12],label=r'$<vw>$')
plt.plot(re2000[:,1],re2000[:,13],'--',label=r'$\mathcal{P}$')
#plt.xlim(0.,1e3)
#plt.ylim(-2.,8.)
lg = plt.legend(loc=1)
plt.xlabel(r'$y^+$')
plt.ylabel(r'Reynolds stresses')
plt.savefig('pb2b')

# Re2000
# y/h, y+, U+, u'+, v'+, w'+, -Om_z+, om_x'+, om_y'+, om_z'+, uv'+, uw'+, vw'+, pr'+, ps'+, psto'+, p'
fig = plt.figure(figsize=(12,8))
plt.plot(re2000[:,1],re2000[:,3],label=r'$<uu>$')
plt.plot(re2000[:,1],re2000[:,4],label=r'$<vv>$')
plt.plot(re2000[:,1],re2000[:,5],label=r'$<ww>$')
plt.plot(re2000[:,1],re2000[:,10],label=r'$<uv>$')
plt.plot(re2000[:,1],re2000[:,11],label=r'$<uw>$')
plt.plot(re2000[:,1],re2000[:,12],label=r'$<vw>$')
plt.plot(re2000[:,1],re2000[:,13],'--',label=r'$\mathcal{P}$')
plt.xlim(0.,400)
#plt.ylim(-2.,8.)
lg = plt.legend(loc=1)
plt.xlabel(r'$y^+$')
plt.ylabel(r'Reynolds stresses')
plt.savefig('pb2b_zoom')
