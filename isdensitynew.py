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
from scipy.integrate import nquad

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

#right_ascension =  #radian(alpha)
declination = np.arange(math.pi * -0.5, math.pi * 0.5, 0.03) #radian(delta)
#distance = 8.4 #kpc

'''
c_icrs = SkyCoord(ra=right_ascension*u.rad, dec=declination*u.rad, distance=distance*u.kpc, frame='icrs')
c_galactic = c_icrs.galactic
'''

'''
print(c)
print(c.icrs)
'''

'''
#check test
c_galactic.representation_type = 'cartesian'
c_icrs.representation_type = 'cartesian'
print(c_galactic)
print(c_icrs)

galactic_x = c_galactic.cartesian.x
galactic_y = c_galactic.cartesian.y
galactic_z = c_galactic.cartesian.z

print(galactic_x, galactic_y, galactic_z) #test get value

print(galactic_x + galactic_y + galactic_z) #test calculation

galacticentric_R = ((8.4*u.kpc - galactic_x) ** 2 + galactic_y ** 2) ** 0.5
galacticentric_z = galactic_z

print(galacticentric_R, galacticentric_z)

galactic_x_value = galactic_x / u.kpc
galactic_y_value = galactic_y / u.kpc
galactic_z_value = galactic_z / u.kpc
galacticentric_R_value = ((8.4 - galactic_x_value) ** 2 + galactic_y_value ** 2) ** 0.5
galacticentric_z_value = galactic_z_value

print(galacticentric_R_value, galacticentric_z_value)
'''

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
    return quad(delta_dependence, 0, 50, args=(delta, alpha), limit=50)[0]

def degree_integral(delta):
    angle_alpha = np.arange(0.05, 2 * math.pi, 0.05)
    sum = distance_integral(delta, 0)
    for i in angle_alpha:
        sum = sum + distance_integral(delta, i)
    return sum

x = declination
y = np.array([degree_integral(t) for t in x])
plt.plot(np.sin(x), y * 3.0857 * 10 ** 21 / 126)
plt.axhline(linewidth = 1, color = 'black')
plt.axvline(linewidth = 1, color = 'black')
plt.show()

'''
# lon_0 is central longitude of projection.
# resolution = 'c' means use crude resolution coastlines.
m = Basemap(projection='hammer',lon_0=0,resolution='c')
m.drawcoastlines()
m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))
m.drawmapboundary(fill_color='aqua')
plt.title("Hammer Projection")
plt.show()
'''
