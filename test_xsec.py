import math
import numpy as np
from params import *
import matplotlib.pyplot as plt

energies = np.linspace(2,10,800)

e_e = lambda e: e - DEL_NP
p_e = lambda e: math.sqrt(e_e(e)**2 - M_E*M_E)
e_exp = lambda e: e**(-0.07056+0.02018*math.log(e)-0.001953*(math.log(e))**3)
xsec = lambda e: 1e-43*p_e(e)*e_e(e)*e_exp(e) # cm^2

xsecs = [xsec(e) for e in energies]

plt.plot(energies,xsecs)
plt.show()
