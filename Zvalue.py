import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

#factor constant
Fp = 0.92
FHe = 0.90
FFe = 0.86

#cross section constant
Betap = 0.082
BetaHe = 0.062
BetaFe = 0.026
sigmapp0 = 17.7 * 10 ** -27 #mb
sigmaHep0 = 60.5 * 10 ** -27 #mb
sigmaFep0 = 551 * 10 ** -27 #mb

#neutrino yield constant
ap0 = 15.5
bp0 = 8.0
Ap = 1

aHe0 = 14.8
bHe0 = 7.3
AHe = 4

aFe0 = 20.0
bFe0 = 5.8
AFe = 56

'''
#CR component flux constant
flux_H = 
flux_He = 
flux_Fe = 
'''

def energy_dependence_p(x, energy):
    energy_x = energy / x

    a = ap0 * (1 + 0.073 * np.math.log(energy_x / 10 ** 6) + 0.0070 * (np.math.log(energy_x / 10 ** 6) ** 2))
    b = bp0 * (1 + 0.020 * np.math.log(energy_x / 10 ** 6) + 0.0018 * (np.math.log(energy_x / 10 ** 6) ** 2))

    ax = Ap * x 

    ax_power = ax ** 0.43
    bax_power = -b * ax_power

    fitting = a * ((1 - ax) ** 3) / x * (math.exp(1) ** bax_power) / (1 + (0.1 / x / energy_x) ** 0.5) ** 2

    cross_section_value = sigmapp0 * (energy_x / 1) ** Betap

    flux_value = 1.3 * (energy_x / 1) ** -2.7

    return fitting / x * cross_section_value * flux_value

def x_integral_p(energy):
    return quad(energy_dependence_p, 0, 1, args=(energy), limit=50)[0]

def energy_dependence_He(x, energy):
    energy_x = energy / x

    a = aHe0 * (1 + 0.073 * np.math.log(energy_x / 10 ** 6) + 0.0070 * (np.math.log(energy_x / 10 ** 6) ** 2))
    b = bHe0 * (1 + 0.020 * np.math.log(energy_x / 10 ** 6) + 0.0018 * (np.math.log(energy_x / 10 ** 6) ** 2))

    ax = AHe * x 

    ax_power = ax ** 0.43
    bax_power = -b * ax_power

    fitting = a * ((1 - ax) ** 3) / x * (math.exp(1) ** bax_power) / (1 + (0.1 / x / energy_x) ** 0.5) ** 2

    cross_section_value = sigmaHep0 * (energy_x / 1) ** BetaHe

    flux_value = 0.54 * (energy_x / 1) ** -2.6

    return fitting / x * cross_section_value * flux_value

def x_integral_He(energy):
    return quad(energy_dependence_He, 0, 1, args=(energy), limit=50)[0]

def x_integral(energy):
    return (Fp * x_integral_p(energy) + FHe * x_integral_He(energy)) * energy ** 2.7

def energy_dependence_He_new(x, energy):
    energy_x = energy / x

    a = aHe0 * (1 + 0.073 * np.math.log(energy_x / 10 ** 6) + 0.0070 * (np.math.log(energy_x / 10 ** 6) ** 2))
    b = bHe0 * (1 + 0.020 * np.math.log(energy_x / 10 ** 6) + 0.0018 * (np.math.log(energy_x / 10 ** 6) ** 2))

    ax = AHe * x 

    ax_power = ax ** 0.43
    bax_power = -b * ax_power

    fitting = a * ((1 - ax) ** 3) / x * (math.exp(1) ** bax_power) / (1 + (0.1 / x / energy_x) ** 0.5) ** 2

    cross_section_value = sigmaHep0 * (energy_x / 1) ** BetaHe

    flux_value = 330 * (energy_x / 1) ** -3

    return fitting / x * cross_section_value * flux_value

def x_integral_He_new(energy):
    return quad(energy_dependence_He_new, 0, 1, args=(energy), limit=50)[0] * FHe * energy ** 2.7

def energy_dependence_p_new(x, energy):
    energy_x = energy / x

    a = ap0 * (1 + 0.073 * np.math.log(energy_x / 10 ** 6) + 0.0070 * (np.math.log(energy_x / 10 ** 6) ** 2))
    b = bp0 * (1 + 0.020 * np.math.log(energy_x / 10 ** 6) + 0.0018 * (np.math.log(energy_x / 10 ** 6) ** 2))

    ax = Ap * x 

    ax_power = ax ** 0.43
    bax_power = -b * ax_power

    fitting = a * ((1 - ax) ** 3) / x * (math.exp(1) ** bax_power) / (1 + (0.1 / x / energy_x) ** 0.5) ** 2

    cross_section_value = sigmapp0 * (energy_x / 1) ** Betap

    flux_value = 330 * (energy_x / 1) ** -3

    return fitting / x * cross_section_value * flux_value

def x_integral_p_new(energy):
    return quad(energy_dependence_p_new, 0, 1, args=(energy), limit=50)[0] * Fp * energy ** 2.7

def energy_dependence_Fe_new(x, energy):
    energy_x = energy / x

    a = aFe0 * (1 + 0.073 * np.math.log(energy_x / 10 ** 6) + 0.0070 * (np.math.log(energy_x / 10 ** 6) ** 2))
    b = bFe0 * (1 + 0.020 * np.math.log(energy_x / 10 ** 6) + 0.0018 * (np.math.log(energy_x / 10 ** 6) ** 2))

    ax = AFe * x 

    ax_power = ax ** 0.43
    bax_power = -b * ax_power

    fitting = a * ((1 - ax) ** 3) / x * (math.exp(1) ** bax_power) / (1 + (0.1 / x / energy_x) ** 0.5) ** 2

    cross_section_value = sigmaFep0 * (energy_x / 1) ** BetaFe

    flux_value = 330 * (energy_x / 1) ** -3

    return fitting / x * cross_section_value * flux_value

def x_integral_Fe_new(energy):
    return quad(energy_dependence_Fe_new, 0, 1, args=(energy), limit=50)[0] * -1 * FFe * energy ** 2.7

x = np.logspace(0, 6.7, 10000, endpoint=True)
y = np.array([x_integral(t) for t in x])
plt.plot(x, y * 0.0077 / (1.6726219 * 10 ** -24))
x1 = np.logspace(6.7, 9.5, 10000, endpoint=True)
y1 = np.array([x_integral_He_new(t) for t in x1])
plt.plot(x1, y1 * 0.0077 / (1.6726219 * 10 ** -24))
x2 = np.logspace(6.7, 9.5, 10000, endpoint=True)
y2 = np.array([x_integral_p_new(t) for t in x2])
plt.plot(x2, y2 * 0.0077 / (1.6726219 * 10 ** -24))
x3 = np.logspace(6.7, 9.5, 10000, endpoint=True)
y3 = np.array([x_integral_Fe_new(t) for t in x3])
plt.plot(x3, y3 * 0.0077 / (1.6726219 * 10 ** -24))
plt.axis([1, 10 ** 8, 5 * 10 ** -8, 5 * 10 ** -5])
plt.xscale('log')
plt.yscale('log')
plt.axhline(linewidth = 1, color = 'black')
plt.axvline(linewidth = 1, color = 'black')
plt.show()