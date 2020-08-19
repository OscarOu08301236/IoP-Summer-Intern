import math
import numpy as np
import matplotlib.pyplot as plt
import astropy as astro
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude # Angles
from mpl_toolkits.basemap import Basemap
from scipy.integrate import quad
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show

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

declination = np.arange(math.pi * -0.5, math.pi * 0.5, 0.03) #radian(delta)
rightascension = np.arange(0, math.pi * 2, 0.03) #radian(alpha)

def delta_dependence(distance, delta, alpha):
    c_icrs = SkyCoord(ra=alpha*u.rad, dec=delta*u.rad, distance=distance*u.kpc, frame='icrs')
    c_galactic = c_icrs.galactic
    c_galactic.representation_type = 'cartesian'

    galactic_x = c_galactic.cartesian.x
    galactic_y = c_galactic.cartesian.y
    galactic_z = c_galactic.cartesian.z

    galactic_x_value = galactic_x / u.kpc
    galactic_y_value = galactic_y / u.kpc
    galactic_z_value = galactic_z / u.kpc
    galacticentric_R_value = ((8.4 - galactic_x_value) ** 2 + galactic_y_value ** 2) ** 0.5
    galacticentric_z_value = galactic_z_value

    R = galacticentric_R_value
    z = np.abs(galacticentric_z_value)

    if R < 7.0:
        ht = ht0 * math.exp(((R - R0) / R0))
        hT = hT0 * math.exp(((R - R0) / R0))
        BRdisk = B0 * math.exp((-1 * (R - RO) / RB))
        grdisk = (BRdisk / B0) ** 0.5

        rhot = proton_mass * nt * math.exp(-(7 - RO) / Rn - z * np.log(2) / ht)
        rhoT = proton_mass * nT * math.exp(-(7 - RO) / Rn - z * np.log(2) / hT)

        rhoh = proton_mass * nh * math.exp(-((R * R + z * z) ** 0.5 - RO) / rh)

        grhalo = (Bh / B0) ** 0.5

        grrho = (rhot + rhoT + rhoh) * grdisk

        return grrho
    elif (R >= 7.0) & (R <= 35.00):
        ht = ht0 * math.exp(((R - R0) / R0))
        hT = hT0 * math.exp(((R - R0) / R0))
        BRdisk = B0 * math.exp((-1 * (R - RO) / RB))
        grdisk = (BRdisk / B0) ** 0.5

        rhot = proton_mass * nt * math.exp(-(R - RO) / Rn - z * np.log(2) / ht)
        rhoT = proton_mass * nT * math.exp(-(R - RO) / Rn - z * np.log(2) / hT)

        rhoh = proton_mass * nh * math.exp(-((R * R + z * z) ** 0.5 - RO) / rh)

        grhalo = (Bh / B0) ** 0.5

        grrho = (rhot + rhoT + rhoh) * grdisk

        return grrho
    elif (R > 35.0) & (R <= 60.0):
        rhoh = proton_mass * nh * math.exp(-((R * R + z * z) ** 0.5 - RO) / rh)

        grhalo = (Bh / B0) ** 0.5

        grrho = rhoh * grhalo
        
        return grrho
    else:
        return 0.0

def distance_integral(delta, alpha):
    return quad(delta_dependence, 0, 45, args=(delta, alpha), limit=45)[0] * 3.0857 * 10 ** 21


x = declination
y = rightascension
z = np.array([[distance_integral(t, s) for t in x] for s in y])
z_average = np.average(z)
maxvalue = np.max(z)
minvalue = np.min(z)
print(z_average)
print(maxvalue, minvalue)
im = imshow(z, interpolation='nearest', aspect='auto'
            , extent=[np.sin(math.pi * -0.5), np.sin(math.pi * 0.5), 0, math.pi * 2])
colorbar(im)
show()
'''
plt.plot(x, y, z * 3.0857 * 10 ** 21)
plt.axhline(linewidth = 1, color = 'black')
plt.axvline(linewidth = 1, color = 'black')
plt.show()
'''
