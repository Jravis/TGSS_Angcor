import numpy as np
import matplotlib.pyplot as plt
from scipy import special
from scipy import integrate

H0 = 70.
C = 3e5 
r0 = 8.
Omega_lam = 0.7
Omega_nr = 0.3
zmin = 0.
zmax = 4.0
amp = 40050.4698
cen = 0.00488450
wid = 0.00866364

amp1 = 133.635376
cen1 = -2.47250015
wid1 = 0.62059123


def ratio(red2):
    return (Omega_lam + Omega_nr*(1.+red2)**3.0)**0.5


def dz_dr(red3):
    return C**-1 * H0 * ratio(red3)


def rint(red1):
    return (C*H0**-1)/(ratio(red1))


def gz(red, ind1, eps1):

    r, err2 = integrate.quad(rint, zmin, zmax)
    return dz_dr(red) * (1. + red)**(ind1-3 - eps1) * r**(1.-ind1)


def dn_dz(z1):

    #return (z1**2) * amp * np.exp(-(z1-cen)**2/ wid)  
    return amp1 * np.exp(-(np.log(z1)-cen1)**2 / wid1)


def int_fun(z, ind, eps):

    return gz(z, ind, eps)*dn_dz(z)**2

    
def correaltion(angle, pw, ep):
    
    #norm_func = lambda zz: (zz**2) * amp * np.exp(-(zz-cen)**2/ wid)
    norm_func = lambda zz: amp1 * np.exp(-(np.log(zz)-cen1)**2/ wid1)  
    intg, err = integrate.quad(int_fun, zmin, zmax,args=(pw, ep) )
    norm, err1 = integrate.quad(norm_func,zmin, zmax)
    B = intg/norm
    fac = special.gamma((indx * 0.5) - 0.5) / special.gamma(indx * 0.5)
    Aw = np.sqrt(np.pi) * r0**indx * fac * B
    return Aw * angle ** (1 - indx)

indx = 1.8

epsi = [indx-3, 0]  # This should roughly the case describing the evolution of 2PCF in theories
# where galaxies identified

col = ['r', 'b', 'y']
theta = 10 ** (np.linspace(np.log10(0.1), np.log10(8.), 100)) # theta range in logrithmic bins 
plt.figure(1, figsize=(8, 6))

for i in xrange(len(epsi)):
    w = correaltion(theta, indx, epsi[i])
    plt.plot(theta, w, '-', color=col[i], linewidth=3)
    plt.xscale("log")
    plt.yscale("log")
    plt.ylim(0.0006, 0.03)
plt.show()

