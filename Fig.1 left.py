import numpy as np 
import matplotlib.pyplot as plt 
import math

#global setting 
E = 10 ** 6 #GeV
a0 = [2.7, 2.7 * 14.8 / 15.5, 2.7 * 20.0 / 15.5]
b0 = [7.7, 7.7 * 7.3 / 8.0, 7.7 * 5.8 / 8.0]
A = [1, 4, 56]

# parameter x
x = np.logspace(-11, 0, 10000, endpoint=True)

#p graph
a1 = a0[0] * (1 + 0.073 * np.math.log(E / 10 ** 6) + 0.0070 * (np.math.log(E / 10 ** 6) ** 2))
b1 = b0[0] * (1 + 0.020 * np.math.log(E / 10 ** 6) + 0.0018 * (np.math.log(E / 10 ** 6) ** 2))

ax1 = A[0] * x 

ax_power1 = ax1 ** 0.43
bax_power1 = -b1 * ax_power1

fitting1 = a1 * ((1 - ax1) ** 3) / x * (math.exp(1) ** bax_power1) / (1 + (0.1 / x / E) ** 0.5) ** 2

#He graph
a2 = a0[1] * (1 + 0.073 * np.math.log(E / 10 ** 6) + 0.0070 * (np.math.log(E / 10 ** 6) ** 2))
b2 = b0[1] * (1 + 0.020 * np.math.log(E / 10 ** 6) + 0.0018 * (np.math.log(E / 10 ** 6) ** 2))

ax2 = A[1] * x 

ax_power2 = ax2 ** 0.43
bax_power2 = -b2 * ax_power2

fitting2 = a2 * ((1 - ax2) ** 3) / x * (math.exp(1) ** bax_power2) / (1 + (0.1 / x / E) ** 0.5) ** 2

#Fe graph
a3 = a0[2] * (1 + 0.073 * np.math.log(E / 10 ** 6) + 0.0070 * (np.math.log(E / 10 ** 6) ** 2))
b3 = b0[2] * (1 + 0.020 * np.math.log(E / 10 ** 6) + 0.0018 * (np.math.log(E / 10 ** 6) ** 2))

ax3 = A[2] * x 

ax_power3 = ax3 ** 0.43
bax_power3 = -b3 * ax_power3

fitting3 = a3 * ((1 - ax3) ** 3) / x * (math.exp(1) ** bax_power3) / (1 + (0.1 / x / E) ** 0.5) ** 2

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

