"""
In this routine we Crossmatch TGSS sources with SDSS and 2df surveys
From SDSS We are using three catalouge Quasar, Galaxie
"""
import numpy as np
import matplotlib.pyplot as plt
from astroML.crossmatch import crossmatch_angular
import matplotlib.gridspec as gridspec

# Read values from 2df best observation data

filename_2df = "best.observations.idz"
Ra_HH = np.genfromtxt(filename_2df, usecols=10)
Ra_MM = np.genfromtxt(filename_2df, usecols=11)
Ra_SS = np.genfromtxt(filename_2df, usecols=12)
Dec_dd = np.genfromtxt(filename_2df, usecols=13)
Dec_MM = np.genfromtxt(filename_2df, usecols= 14)
Dec_SS = np.genfromtxt(filename_2df, usecols= 15)
TwoDF_redsft = np.genfromtxt(filename_2df, usecols= 23)

# convert Ra from HH:MM:SS to degree and Dec dd:MM:SS to degree
Ra_TwoDF = 15.*(Ra_HH + Ra_MM/60. + Ra_SS/3600.)
Dec_TwoDF = (Dec_dd + Dec_MM/60. + Dec_SS/3600.)

# read values from SDSS file having spectroscopic redshift and both galaxies and quasar

filename_AllSdss = "dr12_queries.csv"
Ra = np.genfromtxt(filename_AllSdss, usecols=0, delimiter=',')
Dec = np.genfromtxt(filename_AllSdss, usecols=1, delimiter=',')
u = np.genfromtxt(filename_AllSdss, usecols=2, delimiter=',')
g = np.genfromtxt(filename_AllSdss, usecols=3, delimiter=',')
r = np.genfromtxt(filename_AllSdss, usecols=4, delimiter=',')
red = np.genfromtxt(filename_AllSdss, usecols=5, delimiter=',')
photo_red = np.genfromtxt(filename_AllSdss, usecols=7, delimiter=',')


# TGSS data with ra and dec after applying basic selection criteria flux greater than 3.5 or 4 miliJy

fname_TGSS = "/home/sandeep/CUTE-master/Angcor/out_data/data_17April17/CrossMatch_TGSS_SDSS_50mJy.txt"
RA = np.genfromtxt(fname_TGSS, usecols=0, delimiter='\t')
DEC = np.genfromtxt(fname_TGSS, usecols=1, delimiter='\t')
flux = np.genfromtxt(fname_TGSS, usecols=3, delimiter='\t')

# read values from SDSS file having spectroscopic redshift only galaxies
fname_Gal = "dr12_queries_Galaxies.csv"

Ra_G = np.genfromtxt(fname_Gal, usecols=0, delimiter=',')
Dec_G = np.genfromtxt(fname_Gal, usecols=1, delimiter=',')
u_G = np.genfromtxt(fname_Gal, usecols=2, delimiter=',')
g_G = np.genfromtxt(fname_Gal, usecols=3, delimiter=',')
r_G = np.genfromtxt(fname_Gal, usecols=4, delimiter=',')
red_G = np.genfromtxt(fname_Gal, usecols=7, delimiter=',')
photo_red_G = np.genfromtxt(fname_Gal, usecols=9, delimiter=',')

# read values from SDSS file having spectroscopic redshift only galaxies
fname_Q = "dr12_queries_Quasar.csv"

Ra_Q = np.genfromtxt(fname_Q, usecols=0, delimiter=',')
Dec_Q = np.genfromtxt(fname_Q, usecols=1, delimiter=',')
red_Q = np.genfromtxt(fname_Q, usecols=7, delimiter=',')

data_SDSS = np.zeros((len(Ra), 2), dtype=np.float64)
data_SDSS_G = np.zeros((len(Ra_G), 2), dtype=np.float64) 
data_SDSS_Q = np.zeros((len(Ra_Q), 2), dtype=np.float64) 
data_TGSS = np.zeros((len(RA), 2), dtype=np.float64) 
data_TwoDF = np.zeros((len(Ra_TwoDF), 2), dtype=np.float64) 

data_TGSS[:, 0] = RA
data_TGSS[:, 1] = DEC

data_SDSS[:, 0] = Ra
data_SDSS[:, 1] = Dec

data_TwoDF[:, 0] = Ra_TwoDF
data_TwoDF[:, 1] = Dec_TwoDF

data_SDSS_G[:, 0] = Ra_G
data_SDSS_G[:, 1] = Dec_G

data_SDSS_Q[:, 0] = Ra_Q
data_SDSS_Q[:, 1] = Dec_Q

# crossmatch catalogs taken from astroml by Jake VanderPlas

max_radius = 25./3600  # arcsec

dist, ind = crossmatch_angular(data_TGSS, data_SDSS, max_radius)
dist1, ind1 = crossmatch_angular(data_TGSS, data_TwoDF, max_radius)
dist_G, indG = crossmatch_angular(data_TGSS, data_SDSS_G, max_radius)
dist_Q, indQ = crossmatch_angular(data_TGSS, data_SDSS_Q, max_radius)

# crossmatch_angular gives you distance and index of second array given e.g. data_SDSS
# if not matched then dist for that index is np.inf

print len(RA), len(Ra), dist.shape, ind.shape
print len(RA), len(Ra), dist1.shape, ind1.shape

name1 = '/home/sandeep/CUTE-master/Angcor/SDSSDR12/crossmatch_TGSS_SDSS.txt'
name2 = '/home/sandeep/CUTE-master/Angcor/SDSSDR12/crossmatch_TGSS_TwoDF.txt'
name3 = '/home/sandeep/CUTE-master/Angcor/SDSSDR12/crossmatch_TGSS_SDSS_Galaxies.txt'
name4 = '/home/sandeep/CUTE-master/Angcor/SDSSDR12/crossmatch_TGSS_SDSS_Quasar.txt'

zdis = []
U = []
R = []
G = []
Radio_flux = []

with open(name1, 'w') as f:
    for i in xrange(len(RA)):
        if dist[i] != np.inf:
            zdis.append(red[ind[i]])
            f.write("%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\n" % (RA[i],
                    DEC[i], Ra[ind[i]], Dec[ind[i]], red[ind[i]]))

zdis = np.asarray(zdis)

zdis1 = []
with open(name2, 'w') as f:
    for i in xrange(len(RA)):
        if dist1[i] != np.inf:
            zdis1.append(TwoDF_redsft[ind1[i]])
            f.write("%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\n" % (RA[i],
                    DEC[i], Ra_TwoDF[ind1[i]], Dec_TwoDF[ind1[i]], TwoDF_redsft[ind1[i]]))

zdis1 = np.asarray(zdis1)

zdis2 = []
with open(name3, 'w') as f:
    for i in xrange(len(RA)):
        if dist_G[i] != np.inf:
            zdis2.append(red_G[indG[i]])
            U.append(u_G[indG[i]])
            R.append(r_G[indG[i]])
            G.append(g_G[indG[i]])
            Radio_flux.append(flux[i])
            f.write("%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\n" % (RA[i],
                    DEC[i], Ra_G[indG[i]], Dec_G[indG[i]], red_G[indG[i]],
                    u[indG[i]], r[indG[i]], g[indG[i]], flux[i]))

zdis2 = np.asarray(zdis2)
U = np.asarray(U)
R = np.asarray(R)
G = np.asarray(G)
Radio_flux = np.asarray(Radio_flux)

zdis3 = []
with open(name2, 'w') as f:
    for i in xrange(len(RA)):
        if dist_Q[i] != np.inf:
            zdis3.append(red_Q[indQ[i]])
            f.write("%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\n" % (RA[i],
                    DEC[i], Ra_Q[indQ[i]], Dec_Q[indQ[i]], red_Q[indQ[i]]))

zdis3 = np.asarray(zdis3)


match = ~np.isinf(dist1)
dist_match = dist[match]
dist_match *= 3600
index = (zdis1 < 0.00)
zdis1[index] = 0.0

print min(zdis2), max(zdis2)

bin_l = 10**np.linspace(np.log10(0.0001), np.log10(max(zdis)), 100)
bin_G = 10**np.linspace(np.log10(0.0001), np.log10(max(zdis2)), 100)
bin_Q = 10**np.linspace(np.log10(0.0001), np.log10(max(zdis3)), 100)


plt.figure(1, figsize=(8, 8))

gs = gridspec.GridSpec(2, 2)
ax1 = plt.subplot(gs[0, :])
ax1.hist(zdis, bins=bin_l, histtype='stepfilled', ec='k', fc='#AAAAAA', label='SDSS')
ax1.hist(zdis1, bins=bin_l, histtype='stepfilled', ec='k', fc='dodgerblue', label='2DF')
plt.xlabel('Z', fontsize='large', fontstyle='italic', weight='medium')
plt.ylabel('N', fontsize='large', fontstyle='italic', weight='medium')
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=11)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=11)
plt.xscale("log")
plt.legend()

ax2 = plt.subplot(gs[1, :-1])
ax2.hist(zdis2, bins=bin_G, histtype='stepfilled', ec='k', fc='indianred', label='SDSS_G')
plt.xlabel('Z', fontsize='large', fontstyle='italic', weight='medium')
plt.ylabel('N', fontsize='large', fontstyle='italic', weight='medium')
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=11)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=11)
plt.xscale("log")
plt.legend()


ax3 = plt.subplot(gs[-1:, 1])

ax3.hist(zdis3, bins=bin_Q, histtype='stepfilled', ec='k', fc='lightsalmon', label='SDSS_Q')
plt.xlabel('Z', fontsize='large', fontstyle='italic', weight='medium')
plt.ylabel('N', fontsize='large', fontstyle='italic', weight='medium')
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2,
                labelsize=11)

plt.tick_params(axis='both', which='major', length=8, width=2,
                labelsize=11)
plt.xscale("log")
#plt.yscale("log")
plt.legend()
plt.savefig("cross.eps", dpi=100)

# ********************************************************************************

X = G-R
index = X < 1.0
index1 = X >= 1.0

plt.figure(2, figsize=(6,6))   
plt.plot(Radio_flux[index], X[index], "ro",ms=5, label="Late")
plt.plot(Radio_flux[index1], X[index1], "bo",ms=5, label="Early")
plt.axhline(y=1.0,linewidth=2, linestyle='-.')
plt.xscale("log")
#plt.yscale("log")
plt.legend()
plt.xlabel('F(mJy)', fontsize='large', fontstyle='italic', weight='medium')
plt.ylabel('g-r', fontsize='large', fontstyle='italic', weight='medium')
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)
plt.tight_layout()
plt.ylim(-1,)
plt.savefig("Flux_color.eps", dpi=100)
plt.show()
