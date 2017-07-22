"""
This routine is to create data and 1000 random catalog
"""

import numpy as np
import astropy.io.ascii as ascii
import healpy as hp
import random
from multiprocessing import Process

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
    #name = "/dataspace/sandeep/Angcor/TGSS_data/data_Cat/TGSS_%0.1fmJy.txt" % flux_cut

    #np.savetxt(name, zip(ra, dec, z, flux), fmt='%f', delimiter='\t', newline='\n')

    # This is the section for Random Catalog

    if ran_cat:

        nsample = 50*len(ra)
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
            #name = "/dataspace/sandeep/Angcor/TGSS_data/random_Cat/ranCat_CUTE_%d.txt" % nn
            name = "/dataspace/sandeep/Angcor/TGSS_data/ranCat_CUTE_100times_originaldata.txt"

            print name
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
        nmin = 0
        nmax = 0
        count = 0
        min_core = 1
        max_core = 1#20
        increment =1 #50
        str = []
        for i in xrange(min_core, max_core + 1):
            s = 'Cell_Count%d' % i
            str.append(s)
        print len(str)

        for i in xrange(len(str)):
            nmin = count
            nmax = count + increment
            if nmax == 1000:
                nmax = 1001
            print nmin, nmax, i
            str[i] = Process(target=analysis, args=(rancat, nmin, nmax, cut))
            str[i].start()
            count = nmax

        for i in xrange(len(str)):
            str[i].join()
    else:

        rancat = False
        analysis(rancat, 1, 1, cut)

