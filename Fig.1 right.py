import numpy as np 
import matplotlib.pyplot as plt 
import math

#global setting 
E = [10 ** 4, 10 ** 6, 10 ** 8]
a0 = 2.7
b0 = 7.7
A = 1

# parameter x
x = np.logspace(-11, 0, 10000, endpoint=True)

#10 ^ 4 GeV graph
a1 = a0 * (1 + 0.073 * -2 + 0.0070 * 4)
b1 = b0 * (1 + 0.020 * -2 + 0.0018 * 4)

ax1 = A * x 

ax_power1 = ax1 ** 0.43
bax_power1 = -b1 * ax_power1

fitting1 = a1 * ((1 - ax1) ** 3) / x * (math.exp(1) ** bax_power1) / (1 + (0.1 / x / E[0]) ** 0.5) ** 2

#10 ^ 6 GeV graph
a2 = a0 * (1 + 0.073 * np.math.log(E[1] / 10 ** 6) + 0.0070 * (np.math.log(E[1] / 10 ** 6) ** 2))
b2 = b0 * (1 + 0.020 * np.math.log(E[1] / 10 ** 6) + 0.0018 * (np.math.log(E[1] / 10 ** 6) ** 2))

ax2 = A * x 

ax_power2 = ax2 ** 0.43
bax_power2 = -b2 * ax_power2

fitting2 = a2 * ((1 - ax2) ** 3) / x * (math.exp(1) ** bax_power2) / (1 + (0.1 / x / E[1]) ** 0.5) ** 2

#10 ^ 8 GeV graph
a3 = a0 * (1 + 0.073 * np.math.log(E[2] / 10 ** 6) + 0.0070 * (np.math.log(E[2] / 10 ** 6) ** 2))
b3 = b0 * (1 + 0.020 * np.math.log(E[2] / 10 ** 6) + 0.0018 * (np.math.log(E[2] / 10 ** 6) ** 2))

ax3 = A * x 

ax_power3 = ax3 ** 0.43
bax_power3 = -b3 * ax_power3

fitting3 = a3 * ((1 - ax3) ** 3) / x * (math.exp(1) ** bax_power3) / (1 + (0.1 / x / E[2]) ** 0.5) ** 2

#plotting code
plt.plot(x, x * fitting1)
plt.plot(x, x * fitting2)
plt.plot(x, x * fitting3)
plt.axis([3.7 * 10 ** -11, 1, 0.001, 10])
plt.xscale('log')
plt.yscale('log')
plt.axhline(linewidth = 1, color = 'black')
plt.axvline(linewidth = 1, color = 'black')
plt.show()

