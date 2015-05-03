import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

bur = np.loadtxt('Burgers.out')

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(bur[:,0],bur[:,1])
ax.plot(bur[:,0],bur[:,2])



