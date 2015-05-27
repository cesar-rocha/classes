import numpy as np
import matplotlib.pyplot as plt

from netCDF4 import Dataset

data = Dataset("model_output.nc")


t = data.variables['time'][:]
ke = data.variables['ke'][:]
ens = data.variables['ens'][:]
q = data.variables['q'][:]

