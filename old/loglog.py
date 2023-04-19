# %% codecell

import numpy as np
import matplotlib.pyplot as plt

xx = np.logspace(-2,2,1000)

#plt.loglog(xx,xx)
#plt.loglog(xx,10*xx)
#plt.loglog(xx,100*xx)


D_dig_li = 0.14383 * 1e-5 #DIG cm^2/s
D_tmp_li = 0.14006 * 1e-5 #TMP cm^2/s

D_dig_tfsi = 0.146286 #DIG *(1e-5 cm^2/s)
D_tmp_tfsi = 0.144115 #TMP *(1e-5 cm^2/s)

D_li = 5 * 1e-7
cm2s_to_A2ps = 1e7

MSD_0 = 6 * D_li * cm2s_to_A2ps * xx
MSD_1 = 6 * 0.5  * xx
MSD_3 = 6 * 20   * xx

plt.loglog(xx,   xx, 'g')
plt.loglog(xx,MSD_0, 'r')
plt.loglog(xx,MSD_1, 'b')
plt.loglog(xx,MSD_3, 'orange')

plt.xlim([5e-2,1e2])
plt.ylim([5e0,5e4])
