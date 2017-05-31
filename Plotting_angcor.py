import matplotlib.pyplot as plt
import numpy as np


flux = 50
esti_bis_1 = np.zeros((480, 49), dtype=np.float64)
name = "/dataspace/sandeep/Angcor/TGSS_data/Corr_data/angCor_%dmJy_%d.txt" % (flux, 1)
theta = np.genfromtxt(name, usecols=0, delimiter=',')
for fn in xrange(1, 481):
    name = "/dataspace/sandeep/Angcor/TGSS_data/Corr_data/angCor_%dmJy_%d.txt" % (flux, fn)
    tmp = np.genfromtxt(name, usecols=1, delimiter=',')
    esti_bis_1[fn-1, :] = tmp


mean = np.mean(esti_bis_1, 0, dtype=np.float64)
std_dev = np.std(esti_bis_1, 0, dtype=np.float64)


plt.figure(3, figsize=(8, 6))
plt.plot(theta, mean, '-', color='orange', linewidth=2, label='mean')
plt.fill_between(theta, (mean - std_dev),  (mean + std_dev), alpha=0.5, edgecolor='c',
                 facecolor='paleturquoise')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.minorticks_on()
plt.xlim(0.1,)
plt.ylim(0.001, 0.05)
plt.grid(which='Both')
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=10)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=10)
plt.xlabel(r"$\theta$", fontsize=18)
plt.ylabel(r"$\omega(\theta)$", fontsize=18)
plt.show()
