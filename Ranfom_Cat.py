"""
This routine is to create data and 1000 random catalog
"""

import numpy as np
import astropy.io.ascii as ascii
import healpy as hp
import random
#from numpy import random, pi
from multiprocessing import Process
from numba import jit as njit


seed = 55555333
random.seed(seed)


def rand_sphere(n):
    z = np.arccos(2 * random.rand(n) - 1) - np.pi/2.0  # uniform in -1, 1
    t = 2 * np.pi * random.rand(n)   # uniform in 0, 2*pi
    return t, z


def analysis(ran_cat, nmin, nmax, flux_cut):

    fname_TGSS = "/dataspace/sandeep/Angcor/input/TGSSADR1_7sigma_catalog.tsv"

    data = ascii.read(fname_TGSS, delimiter='\t')

    Ra = data['RA']
    Dec = data['DEC']
    PF = data['Peak_flux']
    Noise = data['RMS_noise']
    Total_flux = data['Total_flux']
    SorceCode = data['Source_code']

    phi, theta = Ra, Dec
    theta = np.abs(np.subtract(theta, 90.))
    NSIDE = 1024  # Since TGSS have 25 arcmin resolution

    print len(Noise), hp.nside2npix(NSIDE)
    print""

    npix = hp.nside2npix(NSIDE)
    healpix_map = np.zeros(hp.nside2npix(NSIDE), dtype=np.double)
    Source_map = np.zeros(hp.nside2npix(NSIDE), dtype=np.double)
    Random_map = np.zeros(hp.nside2npix(NSIDE), dtype=np.double)

    ra = []
    dec = []
    flux = []
    for i in xrange(len(theta)):
        pix = hp.pixelfunc.ang2pix(NSIDE, np.radians(theta[i]), np.radians(phi[i]))
        if Noise[i] < 4.0:
            healpix_map[pix] = 1.0
            if PF[i] >= 32:  # and SorceCode[i]=='S':
                if Total_flux[i] > flux_cut:
                    ra.append((Ra[i]))
                    dec.append((Dec[i]))
                    Source_map[pix] += 1
                    flux.append(Total_flux[i])

    z = np.zeros(len(dec))
    print len(ra)
    name = "/dataspace/sandeep/Angcor/TGSS_data/data_Cat/TGSS_50mJy.txt"
    np.savetxt(name, zip(ra, dec, z, flux), fmt='%f', delimiter='\t', newline='\n')

    # This is the section for Random Catalog

    if ran_cat:

        nsample = 5*len(ra)
#        print nsample
        for nn in xrange(nmin, nmax):

            Ra_rand = []
            Dec_rand = []
            Count = 0

            while Count < nsample:

                theta = np.arccos((1 - 2*random.uniform(0, 1)))
                phi = 2.*np.pi*random.uniform(0, 1)
                ipix = hp.pixelfunc.ang2pix(NSIDE, theta, phi, nest=False)

                if healpix_map[ipix] == 1.0:
                    theta = np.pi/2. - theta
                    Ra_rand.append(np.degrees(phi))
                    Dec_rand.append(np.degrees(theta))
                    Random_map[ipix] += 1
                    Count += 1

#            print len(Ra_rand)
            z = np.zeros(len(Ra_rand))
            name = "/dataspace/sandeep/Angcor/TGSS_data/random_Cat/ranCat_CUTE_%d.txt" % nn
            np.savetxt(name, zip(Ra_rand, Dec_rand, z), fmt='%f', delimiter='\t', newline='\n')
            del Ra_rand
            del Dec_rand
            del z

if __name__ == "__main__":

    print "Enter which flux cut you want"
    cut = float(raw_input(""))

    print "Do you want random catalog Y/N "
    key = raw_input("")
    if key == 'Y':

        rancat = True

        Cell_Count1 = Process(target=analysis, args=(rancat, 1, 31, cut))
        Cell_Count1.start()
        Cell_Count2 = Process(target=analysis, args=(rancat, 31, 61, cut))
        Cell_Count2.start()
        Cell_Count3 = Process(target=analysis, args=(rancat, 61, 91, cut))
        Cell_Count3.start()
        Cell_Count4 = Process(target=analysis, args=(rancat, 91, 121, cut))
        Cell_Count4.start()
        Cell_Count5 = Process(target=analysis, args=(rancat,  121, 151, cut))
        Cell_Count5.start()
        Cell_Count6 = Process(target=analysis, args=(rancat, 151, 181, cut))
        Cell_Count6.start()
        Cell_Count7 = Process(target=analysis, args=(rancat, 181, 211, cut))
        Cell_Count7.start()
        Cell_Count8 = Process(target=analysis, args=(rancat, 211, 241, cut))
        Cell_Count8.start()
        Cell_Count9 = Process(target=analysis, args=(rancat, 241, 271, cut))
        Cell_Count9.start()
        Cell_Count10 = Process(target=analysis, args=(rancat, 271, 301, cut))
        Cell_Count10.start()
        Cell_Count11 = Process(target=analysis, args=(rancat, 301, 331, cut))
        Cell_Count11.start()
        Cell_Count12 = Process(target=analysis, args=(rancat, 331, 361, cut))
        Cell_Count12.start()
        Cell_Count13 = Process(target=analysis, args=(rancat, 361, 391, cut))
        Cell_Count13.start()
        Cell_Count14 = Process(target=analysis, args=(rancat, 391, 421, cut))
        Cell_Count14.start()
        Cell_Count15 = Process(target=analysis, args=(rancat, 421, 451, cut))
        Cell_Count15.start()
        Cell_Count16 = Process(target=analysis, args=(rancat, 451, 481, cut))
        Cell_Count16.start()
        Cell_Count17 = Process(target=analysis, args=(rancat, 511, 541, cut))
        Cell_Count17.start()
        Cell_Count18 = Process(target=analysis, args=(rancat, 541, 571, cut))
        Cell_Count18.start()
        Cell_Count19 = Process(target=analysis, args=(rancat, 571, 601, cut))
        Cell_Count19.start()
        Cell_Count20 = Process(target=analysis, args=(rancat, 601, 631, cut))
        Cell_Count20.start()
        Cell_Count21 = Process(target=analysis, args=(rancat, 631, 661, cut))
        Cell_Count21.start()
        Cell_Count22 = Process(target=analysis, args=(rancat, 661, 691, cut))
        Cell_Count22.start()
        Cell_Count23 = Process(target=analysis, args=(rancat, 691, 721, cut))
        Cell_Count23.start()
        Cell_Count24 = Process(target=analysis, args=(rancat, 721, 751, cut))
        Cell_Count24.start()
        Cell_Count25 = Process(target=analysis, args=(rancat, 751, 781, cut))
        Cell_Count25.start()
        Cell_Count26 = Process(target=analysis, args=(rancat, 781, 831, cut))
        Cell_Count26.start()
        Cell_Count27 = Process(target=analysis, args=(rancat, 831, 861, cut))
        Cell_Count27.start()
        Cell_Count28 = Process(target=analysis, args=(rancat, 861, 900, cut))
        Cell_Count28.start()

        Cell_Count1.join()
        Cell_Count2.join()
        Cell_Count3.join()
        Cell_Count4.join()
        Cell_Count5.join()
        Cell_Count6.join()
        Cell_Count7.join()
        Cell_Count8.join()
        Cell_Count9.join()
        Cell_Count10.join()
        Cell_Count11.join()
        Cell_Count12.join()
        Cell_Count13.join()
        Cell_Count14.join()
        Cell_Count15.join()
        Cell_Count16.join()
        Cell_Count17.join()
        Cell_Count18.join()
        Cell_Count19.join()
        Cell_Count20.join()
        Cell_Count21.join()
        Cell_Count22.join()
        Cell_Count23.join()
        Cell_Count24.join()
        Cell_Count25.join()
        Cell_Count26.join()
        Cell_Count27.join()
        Cell_Count28.join()
    else:

        rancat = False
        analysis(rancat, 1, 1, cut)

