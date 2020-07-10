import numpy as np 
import matplotlib.pyplot as plt 
import math

R = np.arange(-20, 20, 0.02)

B0 = 4 #muG
Bh = 1 #muG
RB = 8.5 #kpc
RO = 8.4 #kpc

BRdisk = B0 * math.exp(1) ** (-1 * (R - RO) / RB)
BRhalo = Bh * math.exp(1) ** (-1 * (R - RO) / RB)
grdisk = (BRdisk / B0) ** 0.5
grhalo = (BRhalo / B0) ** 0.5

plt.plot(R, grdisk)
plt.plot(R, grhalo)
plt.axhline(linewidth = 1, color = 'black')
plt.axvline(linewidth = 1, color = 'black')
plt.show()