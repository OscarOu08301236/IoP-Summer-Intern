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

right_ascension = 0 #radian(alpha)
declination = 0 #radian(delta)
distance = 8.4#kpc

c_icrs = SkyCoord(ra=right_ascension*u.rad, dec=declination*u.rad, distance=distance*u.kpc, frame='icrs')
c_galactic = c_icrs.galactic
print(c_galactic, c_icrs)

c_galactic.representation_type = 'cartesian'
print(c_galactic)

galactic_x = c_galactic.cartesian.x
galactic_y = c_galactic.cartesian.y
galactic_z = c_galactic.cartesian.z
print(galactic_x, galactic_y, galactic_z)

galactic_x_value = galactic_x / u.kpc
galactic_y_value = galactic_y / u.kpc
galactic_z_value = galactic_z / u.kpc
galacticentric_R_value = ((8.4 - galactic_x_value) ** 2 + galactic_y_value ** 2) ** 0.5
galacticentric_z_value = galactic_z_value
print(galactic_x_value, galactic_y_value, galactic_z_value, galacticentric_R_value, galacticentric_z_value)

R = galacticentric_R_value
z = galacticentric_z_value
print(R, z)
