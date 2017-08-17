import numpy as np
import matplotlib.pyplot as plt
import lmfit
import emcee
import corner
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
    B = intg / norm
    fac = special.gamma((pw * 0.5) - 0.5) / special.gamma(pw * 0.5)
    Aw = np.sqrt(np.pi) * R0 ** pw * fac * B
    return Aw * angle ** (1 - pw)

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


dN_dz = []
for i in xrange(len(dn_dz)-1):
    dN_dz.append((dn_dz[i+1]-dn_dz[i])/np.diff(dist_z)[0])


del dn_dz
dn_dz = dN_dz

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

name = "/dataspace/sandeep/Angcor/TGSS_data/Corr_data/50mJy_data/angCor_%dmJy_%d.txt" % (50, 1)
theta = np.genfromtxt(name, usecols=0, delimiter=',')

flx = [200]#, 60, 100, 200]
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


X = []
Y = []

y = mean1[0, :]
x = theta
y_err = std_dev[0, :]

print("Statistics for theta>1.0")
index = (x <= 1.0)  #* (x <= 1.0)
data = np.zeros((3, len(x[index])), dtype=np.double)
data[0, :] = x[index]
data[1, :] = y[index]
data[2, :] = y_err[index]


pars = lmfit.Parameters()
pars.add('r0', value=5.5, min=2.0, max=20.0)
pars.add('indx', value=1.72, min=0.5, max=3.0)

theta1 = x[index]


def residual(p1):
    r = p1['r0'].value
    indx = p1['indx'].value

    epsi = indx - 3  # This should roughly the case describing the evolution of 2PCF in theories
    return data[1, :] - correaltion(theta1, indx, epsi, dn_dz, dist_z, zmin, zmax, r)


def lnprob(p):
    resid = residual(p)
    s = data[2, :]
    resid *= 1 / s
    resid *= resid
    resid += np.log(2 * np.pi * s ** 2)
    return -0.5 * np.sum(resid)

mi = lmfit.Minimizer(residual, pars)
# first solve with Nelder-Mead
#out1 = mi.minimize(method='Nelder')

out1 = mi.minimize(method='diffrential evolution')
mini = lmfit.Minimizer(lnprob, out1.params)

#res = mini.emcee(burn=300, steps=2000, thin=10, nwalkers=200, workers=20, is_weighted=False,
#                 params=mi.params, seed=1230127)

res = mini.emcee(burn=300, steps=2000, thin=10, nwalkers=50, workers=10, is_weighted=False,
                 params=mi.params, seed=1230127)


print('------------------------------------------')
print("median of posterior probability distribution")
print('------------------------------------------')

lmfit.report_fit(res.params)


fig1 = corner.corner(res.flatchain, bins=200, labels=res.var_names,
                    truths=list(res.params.valuesdict().values()), title_fmt=".4f",
                    quantiles=[0.16, 0.5, 0.84], show_titles=True, labels_args={"fontsize": 40})

fig1.set_size_inches(8, 6)
#fig1.savefig('/dataspace/sandeep/Angcor/TGSS_data/limber_data/Photo_Galaxy/PhotoGalaxyemcee_TGSS_fit_200mJy.png', dpi=600)
#fig1.savefig('/dataspace/sandeep/Angcor/TGSS_data/limber_data/Galaxy/Galaxyemcee_TGSS_fit_200mJy.png', dpi=600)


highest_prob = np.argmax(res.lnprob)
hp_loc = np.unravel_index(highest_prob, res.lnprob.shape)
mle_soln = res.chain[hp_loc]
for i, par in enumerate(pars):
    pars[par].value = mle_soln[i]

print("\nMaximum likelihood Estimation")
print('-----------------------------')
print(pars)


plt.figure(2, figsize=(8, 6))
plt.plot(res.flatchain['r0'], 'k-', alpha=0.3)

plt.figure(3, figsize=(8, 6))
plt.plot(res.flatchain['indx'], 'r-', alpha=0.3)





plt.show()