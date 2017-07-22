"""
In this routine we try to fit dn/dz form using parametric way using
Gaussian and lognormal distribution also using non parametric way using
Gaussian regression. We also look for double entry in crossmatch analysis
and rectify it.
"""

import numpy as np
import matplotlib.pyplot as plt
from astroML.crossmatch import crossmatch_angular
from astroML.plotting import hist


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
fname_TGSS = "/dataspace/sandeep/Angcor/TGSS_data/data_Cat/TGSS_50mJy.txt"
RA = np.genfromtxt(fname_TGSS, usecols=0, delimiter='\t')
DEC = np.genfromtxt(fname_TGSS, usecols=1, delimiter='\t')
flux = np.genfromtxt(fname_TGSS, usecols=3, delimiter='\t')

data_TGSS = np.zeros((len(RA), 2), dtype=np.float64)

data_TGSS[:, 0] = RA
data_TGSS[:, 1] = DEC
# crossmatch catalogs
max_radius = 25./3600  # arcsec


print "Enter the key "
key = raw_input("")

if key == 'G':
    fname_Gal = "/dataspace/sandeep/Angcor/SDSSDR12/sdss_data/dr12_queries_Galaxies.csv"
    Ra_G = np.genfromtxt(fname_Gal, usecols=0, delimiter=',')
    Dec_G = np.genfromtxt(fname_Gal, usecols=1, delimiter=',')
    u_G = np.genfromtxt(fname_Gal, usecols=2, delimiter=',')
    g_G = np.genfromtxt(fname_Gal, usecols=3, delimiter=',')
    r_G = np.genfromtxt(fname_Gal, usecols=4, delimiter=',')
    i_G = np.genfromtxt(fname_Gal, usecols=5, delimiter=',')
    z_G = np.genfromtxt(fname_Gal, usecols=6, delimiter=',')
    red_G = np.genfromtxt(fname_Gal, usecols=7, delimiter=',')
    photo_red_G = np.genfromtxt(fname_Gal, usecols=9, delimiter=',')
    data_SDSS = np.zeros((len(Ra_G), 2), dtype=np.float64)
    data_SDSS[:, 0] = Ra_G
    data_SDSS[:, 1] = Dec_G

    dist, ind = crossmatch_angular(data_TGSS, data_SDSS, max_radius)

    zdis = []
    ra_tgss = []
    dec_tgss = []
    ra_sdss = []
    dec_sdss = []
    U = []
    R = []
    G = []
    I = []
    Z = []
    FLUX = []

    for i in xrange(len(RA)):
        if dist[i] != np.inf:
            ra_tgss.append(RA[i])
            dec_tgss.append(DEC[i])
            FLUX.append(flux[i])
            ra_sdss.append(Ra_G[ind[i]])
            dec_sdss.append(Dec_G[ind[i]])
            zdis.append(red_G[ind[i]])
            U.append(u_G[ind[i]])
            R.append(r_G[ind[i]])
            G.append(g_G[ind[i]])
            I.append(i_G[ind[i]])
            Z.append(z_G[ind[i]])

    del data_SDSS
elif key == 'P':

    fname_SDSS_photo = '/dataspace/sandeep/Angcor/SDSSDR12/sdss_data/Photo_Galaxies_sandy_7690.csv'
    Ra_G = np.genfromtxt(fname_SDSS_photo, usecols=1, delimiter=',')
    Dec_G = np.genfromtxt(fname_SDSS_photo, usecols=2, delimiter=',')
    u_G = np.genfromtxt(fname_SDSS_photo, usecols=3, delimiter=',')
    g_G = np.genfromtxt(fname_SDSS_photo, usecols=4, delimiter=',')
    r_G = np.genfromtxt(fname_SDSS_photo, usecols=5, delimiter=',')
    i_G = np.genfromtxt(fname_SDSS_photo, usecols=6, delimiter=',')
    z_G = np.genfromtxt(fname_SDSS_photo, usecols=7, delimiter=',')
    red_G = np.genfromtxt(fname_SDSS_photo, usecols=8, delimiter=',')
    photo_red_Err = np.genfromtxt(fname_SDSS_photo, usecols=9, delimiter=',')
    data_SDSS = np.zeros((len(Ra_G), 2), dtype=np.float64)
    data_SDSS[:, 0] = Ra_G
    data_SDSS[:, 1] = Dec_G

    dist, ind = crossmatch_angular(data_TGSS, data_SDSS, max_radius)

    zdis = []
    ra_tgss = []
    dec_tgss = []
    ra_sdss = []
    dec_sdss = []
    U = []
    R = []
    G = []
    I = []
    Z = []
    FLUX = []

    for i in xrange(len(RA)):
        if dist[i] != np.inf:
            ra_tgss.append(RA[i])
            dec_tgss.append(DEC[i])
            FLUX.append(flux[i])
            ra_sdss.append(Ra_G[ind[i]])
            dec_sdss.append(Dec_G[ind[i]])
            zdis.append(red_G[ind[i]])
            U.append(u_G[ind[i]])
            R.append(r_G[ind[i]])
            G.append(g_G[ind[i]])
            I.append(i_G[ind[i]])
            Z.append(z_G[ind[i]])

    del data_SDSS

elif key == 'Q':
    fname_Q = "/dataspace/sandeep/Angcor/SDSSDR12/sdss_data/dr12_queries_Quasar.csv"
    Ra_Q = np.genfromtxt(fname_Q, usecols=0, delimiter=',')
    Dec_Q = np.genfromtxt(fname_Q, usecols=1, delimiter=',')
    red_Q = np.genfromtxt(fname_Q, usecols=7, delimiter=',')

    data_SDSS_Q = np.zeros((len(Ra_Q), 2), dtype=np.float64)
    data_SDSS_Q[:, 0] = Ra_Q
    data_SDSS_Q[:, 1] = Dec_Q
    dist, ind = crossmatch_angular(data_TGSS, data_SDSS_Q, max_radius)

    ra_tgss = []
    dec_tgss = []
    ra_sdss = []
    dec_sdss = []
    zdis = []
    FLUX = []
    for i in xrange(len(RA)):
        if dist[i] != np.inf:
            ra_tgss.append(RA[i])
            dec_tgss.append(DEC[i])
            FLUX.append(flux[i])
            ra_sdss.append(Ra_Q[ind[i]])
            dec_sdss.append(Dec_Q[ind[i]])
            zdis.append(red_Q[ind[i]])

    del data_SDSS_Q
# ===========================================================================

del data_TGSS
del flux

tgss_ra = []
tgss_dec = []
sdss_ra = []
sdss_dec = []
spec_red = []
sdss_u = []
sdss_r = []
sdss_g = []
sdss_i = []
sdss_z = []
tgss_flux = []
for i in xrange(len(ra_sdss)):
    if ra_sdss[i] not in sdss_ra:
        if dec_sdss[i] not in sdss_dec:
            tgss_ra.append(ra_tgss[i])
            tgss_dec.append(dec_tgss[i])
            tgss_flux.append(FLUX[i])
            sdss_ra.append(ra_sdss[i])
            sdss_dec.append(dec_sdss[i])
            spec_red.append(zdis[i])
            if key == 'G'or key == 'P':
                sdss_i.append(I[i])
                sdss_z.append(Z[i])

del zdis
del ra_tgss
del dec_tgss
del ra_sdss
del dec_sdss
if key == 'G' or key == 'P':
    del U
    del R
    del G
    del I
    del Z
del FLUX

tgss_ra = np.asarray(tgss_ra)
tgss_dec = np.asarray(tgss_dec)
sdss_ra = np.asarray(sdss_ra)
sdss_dec = np.asarray(sdss_dec)
sdss_i = np.asarray(sdss_i)
sdss_z = np.asarray(sdss_z)
spec_red = np.asarray(spec_red)

if key == 'G':
    name = '/dataspace/sandeep/Angcor/SDSSDR12/data/crossmatch_TGSS_SDSS_Galaxies_unique.txt'
elif key == 'Q':
    name = '/dataspace/sandeep/Angcor/SDSSDR12/data/crossmatch_TGSS_SDSS_Quasar_unique.txt'
elif key == 'P':
    name = '/dataspace/sandeep/Angcor/SDSSDR12/data/crossmatch_TGSS_SDSS_Photo_unique.txt'

np.savetxt(name, zip(tgss_ra, tgss_dec, sdss_ra, sdss_dec, spec_red), fmt="%0.6f,%0.6f,%0.6f,%0.6f,%0.6f",
           delimiter=',')


"""
X = sdss_i-sdss_z
plt.figure(1, figsize=(6, 6))
plt.plot(tgss_flux, X, "ro", ms=5, label="Late")
plt.xscale("log")
#plt.yscale("log")
plt.legend()
plt.xlabel('F(mJy)', fontsize='large', fontstyle='italic', weight='medium')
plt.ylabel('i-z', fontsize='large', fontstyle='italic', weight='medium')
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)
plt.tight_layout()
plt.ylim(-1, 2)
#plt.savefig("Flux_color.eps", dpi=100)
"""


if key == 'P':
    index = (spec_red != -9999)
    spec_red = spec_red[index]
print "minimun redshift in crossmatch %f" % min(spec_red)
print "max redshift in crossmatch %f" % max(spec_red)

plt.figure(1, figsize=(8, 6))
h1 = hist(spec_red, bins='knuth', histtype='stepfilled', ec='k', fc='#AAAAAA')
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)
plt.ylabel(r'$N$', fontsize='large', fontstyle='italic', weight='extra bold')
plt.xlabel(r'$Z$', fontsize='large', fontstyle='italic', weight='extra bold')
name = '/dataspace/sandeep/Angcor/SDSSDR12/redshift_histogram%s.eps' % key
plt.savefig(name, dpi=100)
plt.show()



