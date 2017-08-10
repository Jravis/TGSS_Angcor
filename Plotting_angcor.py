import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

"""
name = '/dataspace/sandeep/Angcor/CUTE_TGSS_data/CUTE_Corr_TGSS/CUTE_Corr_100mJy.dat'
theta_100 = np.genfromtxt(name, usecols=0)
xi_100 = np.genfromtxt(name, usecols=1)
xi_err_100 = np.genfromtxt(name, usecols=2)

name = '/dataspace/sandeep/Angcor/CUTE_TGSS_data/CUTE_Corr_TGSS/CUTE_Corr_50mJy.dat'
theta_50 = np.genfromtxt(name, usecols=0)
xi_50 = np.genfromtxt(name, usecols=1)
xi_err_50 = np.genfromtxt(name, usecols=2)

name = '/dataspace/sandeep/Angcor/CUTE_TGSS_data/CUTE_Corr_TGSS/CUTE_Corr_200mJy.dat'
theta_200 = np.genfromtxt(name, usecols=0)
xi_200 = np.genfromtxt(name, usecols=1)
xi_err_200 = np.genfromtxt(name, usecols=2)

name = '/dataspace/sandeep/Angcor/CUTE_TGSS_data/CUTE_Corr_TGSS/CUTE_Corr_60mJy.dat'
theta_60 = np.genfromtxt(name, usecols=0)
xi_60 = np.genfromtxt(name, usecols=1)
xi_err_60 = np.genfromtxt(name, usecols=2)
"""


# K-Tree results

name = "/dataspace/sandeep/Angcor/TGSS_data/Corr_data/50mJy_data/angCor_%dmJy_%d.txt" % (50, 1)
theta = np.genfromtxt(name, usecols=0, delimiter=',')
flx = [50, 60, 100, 200]
mean1 = np.zeros((4, len(theta)), dtype=np.float64)
std_dev = np.zeros((4, len(theta)), dtype=np.float64)

for flux in xrange(len(flx)):
    esti_bis_1 = np.zeros((899, 49), dtype=np.float64)

    for fn in xrange(1, 900):

        if flx[flux] == 50:
            if 481 <= fn < 511:
                name = "/dataspace/sandeep/Angcor/TGSS_data/Corr_data/%dmJy_data/angCor_%d.txt" % (flx[flux], fn)
            else:
                name = "/dataspace/sandeep/Angcor/TGSS_data/Corr_data/%dmJy_data/angCor_%dmJy_%d.txt" % (flx[flux], flx[flux], fn)
        else:
            name = "/dataspace/sandeep/Angcor/TGSS_data/Corr_data/%dmJy_data/angCor_%d.txt" % (flx[flux], fn)
        tmp = np.genfromtxt(name, usecols=1, delimiter=',')
        esti_bis_1[fn-1, :] = tmp

    mean1[flux, :] = np.mean(esti_bis_1, 0, dtype=np.float64)
    std_dev[flux, :] = np.std(esti_bis_1, 0, dtype=np.float64)


plt.figure(1, figsize=(8, 8))
plt.plot(theta, mean1[0, :], '-', color='orange', linewidth=2, label='50mJy')
plt.fill_between(theta, (mean1[0, :] - std_dev[0, :]),  (mean1[0, :] + std_dev[0, :]), alpha=0.5, edgecolor='c',
                 facecolor='#AAAAAA')
plt.plot(theta, mean1[1, :], '-', color='crimson', linewidth=2, label='60mJy')
plt.fill_between(theta, (mean1[1, :] - std_dev[1, :]),  (mean1[1, :] + std_dev[1, :]), alpha=0.5, edgecolor='c',
                 facecolor='#AAAAAA')

plt.plot(theta, mean1[2, :], '-', color='green', linewidth=2, label='100mJy')
plt.fill_between(theta, (mean1[2, :] - std_dev[2, :]),  (mean1[2, :] + std_dev[2, :]), alpha=0.5, edgecolor='c',
                 facecolor='#AAAAAA')

plt.plot(theta, mean1[3, :], '-', color='blue', linewidth=2, label='200mJy')
plt.fill_between(theta, (mean1[3, :] - std_dev[3, :]),  (mean1[3, :] + std_dev[3, :]), alpha=0.5, edgecolor='c',
                 facecolor='#AAAAAA')

plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.minorticks_on()
plt.xlim(0.1,)
plt.ylim(0.001, 0.08)
plt.grid(which='Both')
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=10)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=10)
plt.xlabel(r"$\theta$", fontsize=18)
plt.ylabel(r"$\omega(\theta)$", fontsize=18)

"""
plt.figure(2, figsize=(8, 8))
plt.plot(theta_100, xi_100, '-', color='orange', linewidth=2, label='100mJy')
#plt.fill_between(theta_100, (xi_100 - xi_err_100),  (xi_100 + xi_err_100), alpha=0.5, edgecolor='c',
#                 facecolor='#AAAAAA')
plt.plot(theta_50, xi_50, '-', color='green', linewidth=2, label='50mJy')
#plt.fill_between(theta_50, (xi_50 - xi_err_50),  (xi_50 + xi_err_50), alpha=0.5, edgecolor='c',
#                 facecolor='#AAAAAA')

plt.plot(theta_60, xi_60, '-', color='crimson', linewidth=2, label='60mJy')
#plt.fill_between(theta_60, (xi_60 - xi_err_60),  (xi_60 + xi_err_60), alpha=0.5, edgecolor='c',
#                 facecolor='#AAAAAA')
plt.plot(theta_200, xi_200, '-', color='blue', linewidth=2, label='200mJy')
#plt.fill_between(theta_200, (xi_200 - xi_err_200),  (xi_200 + xi_err_200), alpha=0.5, edgecolor='c',
#                facecolor='#AAAAAA')

plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.minorticks_on()
#plt.xlim(0.1,)
#plt.ylim(0.001, 0.08)
plt.grid(which='Both')
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=10)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=10)
plt.xlabel(r"$\theta$", fontsize=18)
plt.ylabel(r"$\omega(\theta)$", fontsize=18)
"""

plt.figure(3, figsize=(8, 8))
plt.errorbar(theta, mean1[0, :], yerr=std_dev[0, :], fmt='o', color='orange', label='50mJy')
plt.errorbar(theta, mean1[1, :], yerr=std_dev[1, :], fmt='o', color='crimson', label='60mJy')
plt.errorbar(theta, mean1[2, :], yerr=std_dev[2, :], fmt='o', color='green', label='100mJy')
plt.errorbar(theta, mean1[3, :], yerr=std_dev[3, :], fmt='o', color='blue', label='200mJy')

plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.minorticks_on()
plt.xlim(0.1,)
plt.ylim(0.001, 0.08)
plt.grid(which='Both')
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=10)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=10)
plt.xlabel(r"$\theta$", fontsize=18)
plt.ylabel(r"$\omega(\theta)$", fontsize=18)

#plt.savefig("/dataspace/sandeep/Angcor/TGSS_data/AngCorr.png", dpi=600)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

theta_lim = 10 ** np.linspace(np.log10(0.001), np.log10(10), 50)
w = np.zeros((10, len(theta_lim)), dtype=np.float64)

indx = 1.8
r0 = [6.5, 7., 7.2, 7.5, 8.0]
col = ['r', 'b', 'y', 'c', 'g', 'k', 'orange', 'm', 'crimson','lightgreen']
#theta = 10 ** np.linspace(np.log10(0.001), np.log10(10), 50)
epsi = [indx - 3, 0]
count = 0
for r in r0:
    for i in xrange(len(epsi)):
        fname = '/dataspace/sandeep/Angcor/TGSS_data/limber_data/Photo_Galaxy/index_1.8/limber_r0-%0.1f_epsi.%0.1f.txt' % (r, epsi[i])
        w[count, :] = np.genfromtxt(fname, usecols=0, delimiter='\t')
        count+=1


fig = plt.figure(4, figsize=(12, 12))
gs = gridspec.GridSpec(2, 3)

ax1 = plt.subplot(gs[0, 0])
ax1.errorbar(theta, mean1[0, :], yerr=std_dev[0, :], fmt='o', color='orange', label='50mJy')
ax1.errorbar(theta, mean1[1, :], yerr=std_dev[1, :], fmt='o', color='crimson', label='60mJy')
ax1.errorbar(theta, mean1[2, :], yerr=std_dev[2, :], fmt='o', color='green', label='100mJy')
ax1.errorbar(theta, mean1[3, :], yerr=std_dev[3, :], fmt='o', color='blue', label='200mJy')
ax1.plot(theta_lim, w[0, :], '-', color='k', linewidth=2, label=r'($r_{0}$=6.5Mpc, $\epsilon=-1.2$)')
ax1.plot(theta_lim, w[1, :], '-', color='r', linewidth=2, label=r'($r_{0}$=6.5Mpc, $\epsilon=0$)')
ax1.set_xscale('log')
ax1.set_yscale('log')
plt.minorticks_on()
plt.legend(loc=3, prop={'size': 9})
ax1.set_xlim(0.1,)
ax1.set_ylim(0.0001, 0.1)
ax1.set_xlabel(r"$\theta$", fontsize=18)
ax1.set_ylabel(r"$\omega(\theta)$", fontsize=18)


ax2 = plt.subplot(gs[0, 1], sharex=ax1)
ax2.errorbar(theta, mean1[0, :], yerr=std_dev[0, :], fmt='o', color='orange', label='50mJy')
ax2.errorbar(theta, mean1[1, :], yerr=std_dev[1, :], fmt='o', color='crimson', label='60mJy')
ax2.errorbar(theta, mean1[2, :], yerr=std_dev[2, :], fmt='o', color='green', label='100mJy')
ax2.errorbar(theta, mean1[3, :], yerr=std_dev[3, :], fmt='o', color='blue', label='200mJy')
ax2.plot(theta_lim, w[2, :], '-', color='k', linewidth=2, label=r'($r_{0}$=7.0Mpc, $\epsilon=-1.2$)')
ax2.plot(theta_lim, w[3, :], '-', color='r', linewidth=2, label=r'($r_{0}$=7.0Mpc, $\epsilon=0$)')
ax2.set_xscale('log')
ax2.set_yscale('log')
plt.minorticks_on()
ax2.set_xlim(0.1,)
ax2.set_ylim(0.0001, 0.1)
ax2.set_xlabel(r"$\theta$", fontsize=18)
ax2.set_ylabel(r"$\omega(\theta)$", fontsize=18)
plt.legend(loc=3, prop={'size': 9})



ax3 = plt.subplot(gs[0, 2], sharex=ax1)
ax3.errorbar(theta, mean1[0, :], yerr=std_dev[0, :], fmt='o', color='orange', label='50mJy')
ax3.errorbar(theta, mean1[1, :], yerr=std_dev[1, :], fmt='o', color='crimson', label='60mJy')
ax3.errorbar(theta, mean1[2, :], yerr=std_dev[2, :], fmt='o', color='green', label='100mJy')
ax3.errorbar(theta, mean1[3, :], yerr=std_dev[3, :], fmt='o', color='blue', label='200mJy')
ax3.plot(theta_lim, w[4, :], '-', color='k', linewidth=2, label=r'($r_{0}$=7.2Mpc, $\epsilon=-1.2$)')
ax3.plot(theta_lim, w[5, :], '-', color='r', linewidth=2, label=r'($r_{0}$=7.2Mpc, $\epsilon=0 $)')
ax3.set_xscale('log')
ax3.set_yscale('log')
plt.minorticks_on()
ax3.set_xlim(0.1,)
ax3.set_ylim(0.0001, 0.1)
ax3.set_xlabel(r"$\theta$", fontsize=18)
ax3.set_ylabel(r"$\omega(\theta)$", fontsize=18)
plt.legend(loc=3, prop={'size': 9})

ax4 = plt.subplot(gs[1, 0], sharex=ax1)
ax4.errorbar(theta, mean1[0, :], yerr=std_dev[0, :], fmt='o', color='orange', label='50mJy')
ax4.errorbar(theta, mean1[1, :], yerr=std_dev[1, :], fmt='o', color='crimson', label='60mJy')
ax4.errorbar(theta, mean1[2, :], yerr=std_dev[2, :], fmt='o', color='green', label='100mJy')
ax4.errorbar(theta, mean1[3, :], yerr=std_dev[3, :], fmt='o', color='blue', label='200mJy')
ax4.plot(theta_lim, w[6, :], '-', color='k', linewidth=2, label=r'($r_{0}$=7.5Mpc, $\epsilon=-1.2$)')
ax4.plot(theta_lim, w[7, :], '-', color='r', linewidth=2, label=r'($r_{0}$=7.5Mpc, $\epsilon=0 $)')
ax4.set_xscale('log')
ax4.set_yscale('log')
plt.minorticks_on()
ax4.set_xlim(0.1,)
ax4.set_ylim(0.001, 0.08)
ax4.set_xlabel(r"$\theta$", fontsize=18)
ax4.set_ylabel(r"$\omega(\theta)$", fontsize=18)
plt.legend(loc=3, prop={'size': 9})


ax4 = plt.subplot(gs[1, 1], sharex=ax1)
ax4.errorbar(theta, mean1[0, :], yerr=std_dev[0, :], fmt='o', color='orange', label='50mJy')
ax4.errorbar(theta, mean1[1, :], yerr=std_dev[1, :], fmt='o', color='crimson', label='60mJy')
ax4.errorbar(theta, mean1[2, :], yerr=std_dev[2, :], fmt='o', color='green', label='100mJy')
ax4.errorbar(theta, mean1[3, :], yerr=std_dev[3, :], fmt='o', color='blue', label='200mJy')
ax4.plot(theta_lim, w[8, :], '-', color='k', linewidth=2, label=r'($r_{0}$=8.0Mpc, $\epsilon=-1.2$)')
ax4.plot(theta_lim, w[9, :], '-', color='r', linewidth=2, label=r'($r_{0}$=8.0Mpc, $\epsilon=0 $)')
ax4.set_xscale('log')
ax4.set_yscale('log')
plt.minorticks_on()
ax4.set_xlim(0.1,)
ax4.set_ylim(0.001, 0.08)
ax4.set_xlabel(r"$\theta$", fontsize=18)
ax4.set_ylabel(r"$\omega(\theta)$", fontsize=18)
plt.legend(loc=3, prop={'size': 9})


fig.tight_layout() # Or equivalently,  "plt.tight_layout()"
fig.savefig("/dataspace/sandeep/Angcor/TGSS_data/limber_data/Photo_Galaxy/index_1.8/limber.png", dpi=600)


plt.show()
