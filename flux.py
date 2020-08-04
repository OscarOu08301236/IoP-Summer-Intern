import numpy as np 
import matplotlib.pyplot as plt 
import math

energy = np.logspace(0, 8, 10000, endpoint=True)

def energy_result_1(E):
    if (E <= 10 ** 6.7):
        return E ** 2.7 * (3.7 * 10 ** -6 * (E) ** -2.617 + 0.9 * 10 ** -6 * (E) ** -2.538)

    else:
        return E ** 2.7 * (4.4 * 10 ** -4 * (E) ** -2.918)

def energy_result_2(E):
    if (E <= 10 ** 6.7):
        return E ** 2.7 * (3.7 * 10 ** -6 * (E) ** -2.617 + 0.9 * 10 ** -6 * (E) ** -2.538)

    else:
        return E ** 2.7 * (1.2 * 10 ** -4 * (E) ** -2.938)

def energy_result_3(E):
    if (E <= 10 ** 6.7):
        return E ** 2.7 * (3.7 * 10 ** -6 * (E) ** -2.617 + 0.9 * 10 ** -6 * (E) ** -2.538)

    else:
        return E ** 2.7 * (1.3 * 10 ** -5 * (E) ** -2.974)

result_1 = np.array([energy_result_1(t) for t in energy])
result_2 = np.array([energy_result_2(t) for t in energy])
result_3 = np.array([energy_result_3(t) for t in energy])

plt.plot(energy, result_1)
plt.plot(energy, result_2)
plt.plot(energy, result_3)

plt.axis([1, 10 ** 8, 5 * 10 ** -8, 5 * 10 ** -5])
plt.xscale('log')
plt.yscale('log')
plt.axhline(linewidth = 1, color = 'black')
plt.axvline(linewidth = 1, color = 'black')
plt.show()