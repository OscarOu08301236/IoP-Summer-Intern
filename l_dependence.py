import numpy as np 
import matplotlib.pyplot as plt 
import math
from scipy.integrate import quad

#longitude = 0.1

Rn = 3.15 #kpc -> cm
proton_mass = 1.67 * (10 ** -27) * 1000 #g
nt = 1.5 #cm^-3
nT = 0.15 #cm^-3
RO = 8.4 #kpc -> cm
ht0 = 0.15 #kpc -> cm
hT0 = 0.4 #kpc -> cm
R0 = 9.8 #kpc -> cm
nh = 0.0001 #cm^-3
rh = 12 #kpc -> cm
B0 = 4 #muG
Bh = 1 #muG
RB = 8.5  #kpc -> cm

def longitude_dependence(k, longitude):
    R = (8.4 ** 2 + k ** 2 - 2 * 8.4 * k * math.cos(longitude)) ** 0.5

    if R < 7.0:
        BRdisk = B0 * math.exp(1) ** (-1 * (R - RO) / RB)
        grdisk = (BRdisk / B0) ** 0.5

        rhot = proton_mass * nt * math.exp(-(7 - RO) / Rn)
        rhoT = proton_mass * nT * math.exp(-(7 - RO) / Rn)

        rhoh = proton_mass * nh * math.exp(-((R * R) ** 0.5 - RO) / rh)

        grhalo = (Bh / B0) ** 0.5

        grrho = (rhot + rhoT + rhoh) * grdisk

        return grrho
    elif (R >= 7.0) & (R <= 35.00):
        BRdisk = B0 * math.exp(1) ** (-1 * (R - RO) / RB)
        grdisk = (BRdisk / B0) ** 0.5

        rhot = proton_mass * nt * math.exp(-(R - RO) / Rn)
        rhoT = proton_mass * nT * math.exp(-(R - RO) / Rn)

        rhoh = proton_mass * nh * math.exp(-((R * R) ** 0.5 - RO) / rh)

        grhalo = (Bh / B0) ** 0.5

        grrho = (rhot + rhoT + rhoh) * grdisk

        return grrho
    elif (R > 35.0) & (R <= 60.0):
        rhoh = proton_mass * nh * math.exp(-((R * R) ** 0.5 - RO) / rh)

        grhalo = (Bh / B0) ** 0.5

        grrho = rhoh * grhalo
        
        return grrho
    else:
        return 0.0

angle = np.arange(math.pi * -1, math.pi * 1, 0.01)

def k_integral(longitude):
    return quad(longitude_dependence, 0, 50, args=(longitude), limit=50)[0]

#x = np.arange(0, 200, 0.1)

x = angle
y = np.array([k_integral(t) for t in x])
plt.plot(x, y * 3.0857 * 10 ** 21)
plt.axis([-1 * math.pi, math.pi, 0.00, 0.34])
plt.axhline(linewidth = 1, color = 'black')
plt.axvline(linewidth = 1, color = 'black')
plt.show()
