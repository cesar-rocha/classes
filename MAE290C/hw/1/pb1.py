import numpy as np
from numpy import pi
import matplotlib.pyplot as plt


dt = .5
t = np.linspace(0.,10,100)
lam = pi/(2.*dt) + 1j/(2.*dt)

u = np.exp(-lam*t)

fout = np.loadtxt('aoutput.txt')



