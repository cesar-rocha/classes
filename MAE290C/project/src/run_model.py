import numpy as np
import matplotlib.pyplot as plt
import spectral_model as model

reload(model)

m = model.BTModel(nx=128, dt = 0.0001, use_fftw=True)

#m.set_q(np.cos(5.*m.x + 0.*m.y))

m.run()
