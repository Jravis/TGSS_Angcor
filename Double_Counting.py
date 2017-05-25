"""
In this routine we try to fit dn/dz form using parametric way using
Gaussian and lognormal distribution also using non parametric way using
Gaussian regression. We also look for double entry in crossmatch analysis
and rectify it.
"""

import numpy as np
import matplotlib.pyplot as plt
from astroML.crossmatch import crossmatch_angular
from lmfit import Model
from sklearn.gaussian_process import GaussianProcessRegressor, kernels

H0 = 70.
C = 3e5 
r0 = 8.
Omega_lam = 0.7
Omega_nr = 0.3
zmin = 0.
zmax = 4.0


# gaussian model with z**2 variation


def gal_dist(x, amp, cen, wid):
    return (x**2) * amp * np.exp(-(x-cen)**2/ wid)  

# lognormal distribution


def gal_dist1(x, amp1, cen1, wid1):
    return amp1 * np.exp(-(np.log(x)-cen1)**2 / wid1)

# Histogram routine


def Hist(data, mbin):
    frq = []
    for j in xrange(len(mbin)):
        count2 = 0
        for ii in xrange(len(data)):
            if j != (len(mbin)-1):
                if mbin[j] <= data[ii] < mbin[j+1]:
                    count2 += 1
        frq.append(count2)
    frq = np.asarray(frq, dtype=float)
    return frq

fname_SDSS_photo = 'Photo_Galaxies_sandy_7690.csv'
fname_SDSS = "dr12_queries.csv"
Ra = np.genfromtxt(fname_SDSS, usecols=0, delimiter=',')
Dec = np.genfromtxt(fname_SDSS, usecols=1, delimiter=',')
red = np.genfromtxt(fname_SDSS, usecols=5, delimiter=',')

fname_TGSS = "/home/sandeep/CUTE-master/Angcor/out_data/data_17April17/CrossMatch_TGSS_SDSS_50mJy.txt"
RA = np.genfromtxt(fname_TGSS, usecols=0, delimiter='\t')
DEC = np.genfromtxt(fname_TGSS, usecols=1, delimiter='\t')
flux = np.genfromtxt(fname_TGSS, usecols=3, delimiter='\t')

data_SDSS = np.zeros((len(Ra), 2), dtype=np.float64) 
data_TGSS = np.zeros((len(RA), 2), dtype=np.float64) 

data_TGSS[:, 0] = RA
data_TGSS[:, 1] = DEC

data_SDSS[:, 0] = Ra
data_SDSS[:, 1] = Dec

# crossmatch catalogs

max_radius = 25./3600  # arcsec

dist, ind = crossmatch_angular(data_TGSS, data_SDSS, max_radius)
del data_SDSS
del data_TGSS
del flux

name = '/home/sandeep/CUTE-master/Angcor/SDSSDR12/crossmatch_TGSS_SDSS_unique.txt'
zdis = []
ra_tgss = []
dec_tgss = []
ra_sdss = []
dec_sdss = []

for i in xrange(len(RA)):
    if dist[i] != np.inf:
        ra_tgss.append(RA[i])
        dec_tgss.append(DEC[i])
        ra_sdss.append(Ra[ind[i]])
        dec_sdss.append(Dec[ind[i]])
        zdis.append(red[ind[i]])
       
tgss_ra = []
tgss_dec = []
sdss_ra = []
sdss_dec = []
spec_red = []

for i in xrange(len(ra_sdss)):
    if ra_sdss[i] not in sdss_ra:
        if dec_sdss[i] not in sdss_dec:
            tgss_ra.append(ra_tgss[i])
            tgss_dec.append(dec_tgss[i])
            sdss_ra.append(ra_sdss[i])
            sdss_dec.append(dec_sdss[i])
            spec_red.append(zdis[i])

tgss_ra = np.asarray(tgss_ra)
tgss_dec = np.asarray(tgss_dec)
sdss_ra = np.asarray(sdss_ra)
sdss_dec = np.asarray(sdss_dec)
spec_red = np.asarray(spec_red)

np.savetxt(name, zip(tgss_ra, tgss_dec, sdss_ra, sdss_dec, spec_red), fmt="%0.6f,%0.6f,%0.6f,%0.6f,%0.6f",
           delimiter=',')
print "minimun redshift in crossmatch %f" % min(spec_red)
print "max redshift in crossmatch %f" % max(spec_red)

# calculate distribution of redsift from crossmatch data using logarithmic bins in Z
bin_l = 10 ** np.linspace(np.log10(0.0001), np.log10(max(spec_red)), 50)
y = Hist(spec_red, bin_l)

# Using Lmfit we try to fit our two models gal_dist and gal_dist1 to data
# initial guess
wid_g = np.std(bin_l)**2
cen_g = np.mean(bin_l)
amp_g = 100. 

gmodel = Model(gal_dist)
result = gmodel.fit(y, x=bin_l, amp= amp_g, cen = cen_g, wid=wid_g)

gmodel = Model(gal_dist1)
result1 = gmodel.fit(y, x=bin_l, amp1= amp_g, cen1 = cen_g, wid1=wid_g)
 
print result.fit_report()
print result1.fit_report()

# defining gaussian kernal for gaussian regression it's taken form
# David parkinson ipython notebook

kernel = kernels.RBF(length_scale=2.0, length_scale_bounds=(0.1, 10.0))
gp = GaussianProcessRegressor(kernel=kernel)

xprime = np.atleast_2d(bin_l).T
yprime = np.atleast_2d(y).T
gp.fit(xprime, yprime)
xplot_prime = np.atleast_2d(bin_l).T
y_pred, sigma = gp.predict(bin_l[:, np.newaxis], return_std=True)
y_pred = y_pred.flatten()


plt.figure(1, figsize=(8, 8))   
plt.plot(bin_l, y, 'bo')
#plt.plot(bin_l, result.init_fit, 'k--')
plt.plot(bin_l, result.best_fit, 'r-')
plt.plot(bin_l, result1.best_fit, 'g-')
plt.xscale("log")
plt.ylim(-20, 160)
#plt.savefig("model_fit_50bin_log.eps", dpi=100)


plt.figure(2)
h0 = plt.plot(bin_l, y, linestyle='None',marker='s',color='black')
h1 = plt.plot(bin_l, y_pred, label='gaussian process')
plt.fill_between(bin_l, y_pred - sigma, y_pred + sigma,alpha=0.5,
                 color='skyblue')
plt.legend()
plt.ylabel('y')
plt.xlabel('x')
plt.xscale("log")


plt.show()



