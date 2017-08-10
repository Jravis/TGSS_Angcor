"""
We Obtain redshift distrbution of crossmatch catalouge and by fitting non paramatrically
We obtain limber inversion
"""

from scipy import special
from scipy import integrate
import scipy.interpolate
import numpy as np
import matplotlib.pyplot as plt
from astroML.plotting import hist
from scipy.integrate import simps


H0 = 70.
C = 3e5
Omega_lam = 0.7
Omega_nr = 0.3


def ratio(red2):
    return (Omega_lam + Omega_nr * (1. + red2) ** 3.0) ** 0.5


def dz_dr(red3):
    return C**-1 * H0 * ratio(red3)


def rint(red1):
    return (C * H0 ** -1) / (ratio(red1))


def gz(red, ind1, eps1, Zmin, Zmax):
    r, err2 = integrate.quad(rint, Zmin, Zmax)
    return dz_dr(red) * (1. + red) ** (ind1 - 3 - eps1) * r ** (1. - ind1)


def simpson38(f):
    Sum = f[0]+f[len(f)-1]
    for i in xrange(1, len(f)-1):
        if i % 3 == 0:
            Sum += 2*f[i]
        else:
            Sum += 3*f[i]
    return Sum


def int_func(z, ind, eps, distfunc, Zmin, Zmax):

    func = []

    for i in xrange(len(z)-1):
        func.append(gz(z[i], ind, eps, Zmin, Zmax) * distfunc[i]**2)
    h = np.diff(z)
    integral = h[0]*simpson38(func)*3./8.
    #integral = simps(func, z[:-1])
    return integral


def norm_func(z1, disfunc):
    h = np.diff(z1)
    return (h[0]*simpson38(disfunc)*3./8.)**2
#    return (simps(disfunc, z1[:-1]))**2


def correaltion(angle, pw, ep, dN_dZ, zrange , zl, zu, R0):

    intg = int_func(zrange, pw, ep, dN_dZ, zl, zu)
    norm = norm_func(zrange, dN_dZ)
    print intg, norm
    B = intg / norm
    fac = special.gamma((pw * 0.5) - 0.5) / special.gamma(pw * 0.5)
    Aw = np.sqrt(np.pi) * R0 ** pw * fac * B
    return Aw * angle ** (1 - pw)

#  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

print "Enter the key "
key = raw_input("")
if key == 'G':
    fname = '/dataspace/sandeep/Angcor/SDSSDR12/data/crossmatch_TGSS_SDSS_Galaxies_unique.txt'
    redshift = np.genfromtxt(fname, usecols=4, delimiter=',')
elif key == 'Q':
    fname = '/dataspace/sandeep/Angcor/SDSSDR12/data/crossmatch_TGSS_SDSS_Quasar_unique.txt'
    redshift = np.genfromtxt(fname, usecols=4, delimiter=',')
elif key == 'P':
    fname = '/dataspace/sandeep/Angcor/SDSSDR12/data/crossmatch_TGSS_SDSS_Photo_unique.txt'
    redshift = np.genfromtxt(fname, usecols=4, delimiter=',')


index = (redshift >= 0.01)

zmin = min(redshift[index])
zmax = max(redshift[index])


y, bin_edges = np.histogram(redshift[index], bins='auto')

# Using Cubic spline to fit data non parametrically


fint = scipy.interpolate.interp1d(bin_edges[:-1], y, kind='cubic')

dn_dz = []
dist_z = np.arange(min(bin_edges[:-1]), max(bin_edges[:-1]), 0.001)
for i in range(0, len(dist_z)):
    dn_dz.append(fint(dist_z[i]))
dn_dz = np.asarray(dn_dz)


plt.figure(1, figsize=(8, 6))
h1 = hist(redshift[index], bins='knuth', histtype='stepfilled', ec='k', fc='#AAAAAA', label='Knuth')

x_test = h1[1]
x_test = np.asarray(x_test)
fi = scipy.interpolate.interp1d(x_test[:-1], h1[0], kind='cubic')
y_interp =[]
x_interp = np.arange(min(x_test[:-1]), max(x_test[:-1]), 0.001)
for i in range(0, len(x_interp)):
    y_interp.append(fi(x_interp[i]))

plt.plot(x_interp, y_interp, 'r')

plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)
plt.ylabel(r'$N$', fontsize='large', fontstyle='italic', weight='extra bold')
plt.xlabel(r'$Z$', fontsize='large', fontstyle='italic', weight='extra bold')
#plt.xscale("log")
plt.legend(loc=1)
plt.title("Redshift distribution")
#plt.savefig('/dataspace/sandeep/Angcor/TGSS_data/limber_data/QuasarSpecred_knuth_bin_dist.eps', dpi=100)


plt.figure(2, figsize=(8, 6))
#plt.plot(bin_edges[:-1], y, 'bo', label='SpecZ data')
plt.plot(bin_edges[:-1], y, 'bo', label='PhotoZ data')
plt.plot(dist_z, dn_dz, 'r', label='fit')
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)
plt.ylabel(r'$N$', fontsize='large', fontstyle='italic', weight='extra bold')
plt.xlabel(r'$Z$', fontsize='large', fontstyle='italic', weight='extra bold')
plt.xscale("log")
plt.legend()
plt.legend(loc=1)
plt.title("Redshift distribution")
plt.ylim(0, 310)
plt.savefig('/dataspace/sandeep/Angcor/TGSS_data/limber_data/Photo_Galaxy/PhotoGalaxySpecred_autobin_dist.png', dpi=600)


dN_dz = []
for i in xrange(len(dn_dz)-1):
    dN_dz.append((dn_dz[i+1]-dn_dz[i])/np.diff(dist_z)[0])


del dn_dz
dn_dz = dN_dz

indx = 1.8
epsi = [indx - 3, 0]  # This should roughly the case describing the evolution of 2PCF in theories
# where galaxies identified

r0 = [6.5, 7., 7.2, 7.5, 8.0]
col = ['r', 'b', 'y', 'c', 'g', 'k', 'orange', 'm', 'crimson','lightgreen']
theta = 10 ** np.linspace(np.log10(0.001), np.log10(10), 50)
count = 0

plt.figure(3, figsize=(8, 6))
for r in r0:
    for i in xrange(len(epsi)):
        w = correaltion(theta, indx, epsi[i], dn_dz, dist_z, zmin, zmax, r)
        fname = '/dataspace/sandeep/Angcor/TGSS_data/limber_data/Photo_Galaxy/index_1.8/limber_r0-%0.1f_epsi.%0.1f.txt' % (r, epsi[i])
        print fname
        np.savetxt(fname, zip(w, theta), fmt='%0.6e', delimiter='\t')
        plt.plot(theta, w, '-', color=col[count], linewidth=3, label='(%0.1f Mpc, %0.1f)' % (r, epsi[i]))
        plt.xscale("log")
        plt.yscale("log")
    #    plt.ylim(0.0006, 0.03)
        count += 1
plt.legend()
plt.show()
