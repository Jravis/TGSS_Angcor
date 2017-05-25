"""In this code we create random and sample catalog
for TGSS"""

import numpy as np
import matplotlib.pyplot as plt
import math as m
from astroML.crossmatch import crossmatch_angular
from astroML.plotting import hist
import sys
from time import time
import matplotlib.gridspec as gridspec


RA = np.genfromtxt("/home/sandeep/CUTE-master/Angcor/out_data/data_17April17/CrossMatch_TGSS_SDSS_50mJy.txt",
                   usecols= 0, delimiter='\t')
DEC = np.genfromtxt("/home/sandeep/CUTE-master/Angcor/out_data/data_17April17/CrossMatch_TGSS_SDSS_50mJy.txt",
                    usecols= 1, delimiter='\t')
flux = np.genfromtxt("/home/sandeep/CUTE-master/Angcor/out_data/data_17April17/CrossMatch_TGSS_SDSS_50mJy.txt",
                    usecols= 3, delimiter='\t')

Ra_G = np.genfromtxt("Photo_Galaxies_sandy_7690.csv", usecols= 1, delimiter=',')
Dec_G = np.genfromtxt("Photo_Galaxies_sandy_7690.csv", usecols= 2, delimiter=',')
u_G = np.genfromtxt("Photo_Galaxies_sandy_7690.csv", usecols= 3, delimiter=',')
g_G = np.genfromtxt("Photo_Galaxies_sandy_7690.csv", usecols= 4, delimiter=',')
r_G = np.genfromtxt("Photo_Galaxies_sandy_7690.csv", usecols= 5, delimiter=',')
red_G = np.genfromtxt("Photo_Galaxies_sandy_7690.csv", usecols= 6, delimiter=',')
photo_red_Err = np.genfromtxt("Photo_Galaxies_sandy_7690.csv", usecols= 7, delimiter=',')

data_SDSS_G = np.zeros((len(Ra_G), 2), dtype=np.float64) 
data_TGSS = np.zeros((len(RA), 2), dtype=np.float64) 

data_TGSS[:, 0] = RA
data_TGSS[:, 1] = DEC

data_SDSS_G[:, 0] = Ra_G
data_SDSS_G[:, 1] = Dec_G

# crossmatch catalogs

max_radius = 25./3600 #  arcsec
dist_G, indG = crossmatch_angular(data_TGSS, data_SDSS_G, max_radius)
name3 = '/home/sandeep/CUTE-master/Angcor/SDSSDR12/crossmatch_TGSS_Photo_SDSS_Galaxies.txt'

U=[]
R=[]
G=[]
Radio_flux = []

zdis2=[]
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
                    u_G[indG[i]], r_G[indG[i]], g_G[indG[i]], flux[i] ))

zdis2=np.asarray(zdis2)
U=np.asarray(U)
R=np.asarray(R)
G=np.asarray(G)
Radio_flux=np.asarray(Radio_flux)

bin_G = 10 ** np.linspace(np.log10(0.0001), np.log10(max(zdis2)), 100)
plt.figure(1, figsize=(8, 8))   
plt.hist(zdis2, bins=bin_G, histtype='stepfilled', ec='k', fc='indianred', label='SDSS_G')
plt.xlabel('Z', fontsize='large', fontstyle='italic', weight='medium')
plt.ylabel('N', fontsize='large', fontstyle='italic', weight='medium')
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2,
                labelsize=11)
plt.tick_params(axis='both', which='major', length=8, width=2,
                labelsize=11)
plt.xscale("log")
plt.legend()
plt.savefig("Photo_cross.eps", dpi=100)
X = G-R
index = (X)<1.0
index1 = (X)>=1.0

plt.figure(2, figsize=(6,6))   
plt.plot( Radio_flux[index], X[index], "ro",ms=5, label="Late")
plt.plot( Radio_flux[index1], X[index1], "bo",ms=5, label="Early")
plt.axhline(y=1.0,linewidth=2, linestyle='-.')
plt.xscale("log")
#plt.yscale("log")
plt.legend()
plt.xlabel('F(mJy)', fontsize='large', fontstyle='italic', weight='medium')
plt.ylabel('g-r', fontsize='large', fontstyle='italic', weight='medium')
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)
#plt.grid(which='both')
plt.tight_layout()
#plt.ylim(,)
plt.savefig("photo_Flux_color.eps", dpi=100)

plt.show()
