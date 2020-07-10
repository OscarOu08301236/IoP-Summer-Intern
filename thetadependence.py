import numpy as np 
import matplotlib.pyplot as plt 
import math
import astropy as astro
from mpl_toolkits.basemap import Basemap

R7 = np.arange(7, 20, 0.02)
R = np.arange(0, 7, 0.02)
sin_sigma = np.arange(-1.0, 1.0, 0.1)

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

B0 = 4 #muG
Bh = 1 #muG
RB = 8.5 #kpc

BRdisk0 = B0 * math.exp(1) ** (-1 * (R - RO) / RB)
BRhalo0 = Bh * math.exp(1) ** (-1 * (R - RO) / RB)
grdisk0 = (BRdisk0/ B0) ** 0.5
grhalo0 = (BRhalo0 / B0) ** 0.5

BRdisk = B0 * math.exp(1) ** (-1 * (R7 - RO) / RB)
BRhalo = Bh * math.exp(1) ** (-1 * (R7 - RO) / RB)
grdisk = (BRdisk/ B0) ** 0.5
grhalo = (BRhalo / B0) ** 0.5

for i in sin_sigma:
    rhot = proton_mass * nt * np.exp(-(R7 - RO) / Rn - (R7 * i) * np.math.log(2) / ht7)
    rhot0 = proton_mass * nt * np.exp(-(7 - RO) / Rn - (R * i) * np.math.log(2) / ht)
    rhoT = proton_mass * nT * np.exp(-(R7 - RO) / Rn - (R7 * i) * np.math.log(2) / hT7)
    rhoT0 = proton_mass * nT * np.exp(-(7 - RO) / Rn - (R * i) * np.math.log(2) / hT)

    rhoh = proton_mass * nh * np.exp(-((R7 * R7 + (R7 * i) * (R7 * i)) ** 0.5 - RO) / rh)
    rhoh0 = proton_mass * nh * np.exp(-((R * R + (R * i) * (R * i)) ** 0.5 - RO) / rh)

    rhototaldisk = rhot + rhoT 
    rhototal0disk = rhot0 + rhoT0 

    grhodisk = rhototaldisk * grdisk
    grhodisk0 = rhototal0disk * grdisk0

    grhohalo = rhoh * grhalo
    grhohalo0 = rhoh0 * grhalo0

    grho = grhodisk + grhohalo
    grho0 = grhodisk0 + grhohalo0

    plt.plot(R7, grho)
    plt.plot(R, grho0)

plt.axhline(linewidth = 1, color = 'black')
plt.axvline(linewidth = 1, color = 'black')
plt.show()

"""
alpha = -93 #degree
theta = -29 #degree
unit_vector_x = math.cos(theta) * math.cos(alpha)
unit_vector_y = math.cos(theta) * math.sin(alpha)
unit_vector_z = math.sin(theta)
R = np.arange(-20, 7, 0.02)
R7 = np.arange(7, 20, 0.02)
position_x = R * unit_vector_x
position_y = R * unit_vector_y
position_z = R * unit_vector_z
"""