import numpy as np 
import matplotlib.pyplot as plt 
import math

R7 = np.arange(7, 20, 0.02)
R = np.arange(-20, 7, 0.02)
z = np.arange(0.5, 5, 0.5)

Rn = 3.15 #kpc 
proton_mass = 1.67 * (10 ** -27)
nt = 1.5 #cm^-3
nT = 0.15 #cm^-3
RO = 8.4 #kpc
ht0 = 0.15 #kpc
hT0 = 0.4 #kpc
R0 = 9.8 #kpc
ht7 = ht0 * math.exp(1) ** ((R7 - R0) / R0)
ht = ht0 * math.exp(1) ** ((R - R0) / R0)
hT7 = hT0 * math.exp(1) ** ((R7 - R0) / R0)
hT = hT0 * math.exp(1) ** ((R - R0) / R0)
nh = 0.001 #cm^-3
rh = 12 #kpc

for i in z :
    rhot = proton_mass * nt * np.exp(-(R7 - RO) / Rn - i * np.math.log(2) / ht7)
    rhot0 = proton_mass * nt * np.exp(-(7 - RO) / Rn - i * np.math.log(2) / ht)
    rhoT = proton_mass * nT * np.exp(-(R7 - RO) / Rn - i * np.math.log(2) / hT7)
    rhoT0 = proton_mass * nT * np.exp(-(7 - RO) / Rn - i * np.math.log(2) / hT)

    rhoh = proton_mass * nh * np.exp(-((R7 * R7 + i * i) ** 0.5 - RO) / rh)
    rhoh0 = proton_mass * nh * np.exp(-((R * R + i * i) ** 0.5 - RO) / rh)

    rhototal = rhot + rhoT + rhoh
    rhototal0 = rhot0 + rhoT0 + rhoh0

    """
    plt.plot(R7, rhot)
    plt.plot(R, rhot0)
    plt.plot(R7, rhoT)
    plt.plot(R, rhoT0)
    plt.plot(R7, rhoh)
    plt.plot(R, rhoh0)
    """
    plt.plot(R7, rhototal)
    plt.plot(R, rhototal0)

"""
plt.xscale('log')
plt.yscale('log')
"""
plt.axhline(linewidth = 1, color = 'black')
plt.axvline(linewidth = 1, color = 'black')
plt.show()