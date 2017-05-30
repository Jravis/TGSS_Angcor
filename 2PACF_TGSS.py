"""
This routine is to calculate angular correlation function given
Random catalog and uses KD-tree method by Jake VanderPals
"""


import numpy as np
from sklearn.neighbors import BallTree
from multiprocessing import Process
import random

seed = 55555333

random.seed(55555333)

# Check if scikit-learn's two-point functionality is available.
# This was added in scikit-learn version 0.14
try:
    from sklearn.neighbors import KDTree
    sklearn_has_two_point = True
except ImportError:
    import warnings
    sklearn_has_two_point = False


def angular_dist_to_euclidean_dist(D, r=1):

    """convert angular distances to euclidean
    distances"""
    return 2 * r * np.sin(0.5 * D * np.pi / 180.)


def ra_dec_to_xyz(ra, dec):

    """
    We convert (ra, dec) to (x,y,z) on a sphere
    :param ra: of galaxy
    :param dec: of galaxy
    :return: (x,y,z)
    """
    sin_ra = np.sin(ra * np.pi / 180.)
    cos_ra = np.cos(ra * np.pi / 180.)

    sin_dec = np.sin(np.pi / 2 - dec * np.pi / 180.)
    cos_dec = np.cos(np.pi / 2 - dec * np.pi / 180.)

    return (cos_ra * sin_dec,
            sin_ra * sin_dec,
            cos_dec)


def two_point(data, bins, method='standard', data_R=None, random_state=None):

    """
    :param data: (x,y,z) data of galaxies
    :param bins: euclidean distance bins
    :param method: Landy-Szalay estimator
    :param data_R: random data
    :param random_state:
    :return:
    """
    if method not in ['standard', 'landy-szalay']:
        raise ValueError("method must be 'standard' or 'landy-szalay'")
    if bins.ndim != 1:
        raise ValueError("bins must be a 1D array")
    if data.ndim == 1:
        data = data[:, np.newaxis]
    elif data.ndim != 2:
        raise ValueError("data should be 1D or 2D")

    n_samples, n_features = data.shape

    Nbins = len(bins) - 1

    # shuffle all but one axis to get background distribution

    factor = len(data_R) * 1. / len(data)
    if sklearn_has_two_point:
        # Fast two-point correlation functions added in scikit-learn v. 0.14
        KDT_D = KDTree(data)
        KDT_R = KDTree(data_R)
        counts_DD = KDT_D.two_point_correlation(data, bins)
        counts_RR = KDT_R.two_point_correlation(data_R, bins)

    else:
        warnings.warn("Version 0.3 of astroML will require scikit-learn "
                      "version 0.14 or higher for correlation function "
                      "calculations. Upgrade to sklearn 0.14+ now for much "
                      "faster correlation function calculations.")

        BT_D = BallTree(data)
        BT_R = BallTree(data_R)

        counts_DD = np.zeros(Nbins + 1)
        counts_RR = np.zeros(Nbins + 1)

        for i in range(Nbins + 1):
            counts_DD[i] = np.sum(BT_D.query_radius(data, bins[i], count_only=True))
            counts_RR[i] = np.sum(BT_R.query_radius(data_R, bins[i], count_only=True))

    DD1 = np.diff(counts_DD)
    RR1 = np.diff(counts_RR)

    DD = counts_DD
    RR = counts_RR

    RR_zero = (RR == 0)

    RR[RR_zero] = 1

    if sklearn_has_two_point:
        counts_DR = KDT_R.two_point_correlation(data, bins)
    else:
        counts_DR = np.zeros(Nbins + 1)
        for i in range(Nbins + 1):
            counts_DR[i] = np.sum(BT_R.query_radius(data, bins[i], count_only=True))

    DR1 = np.diff(counts_DR)
    DR = counts_DR
    corr = (factor ** 2 * DD - 2 * factor * DR + RR) / RR
    corr1 = (factor ** 2 * DD1 - 2 * factor * DR1 + RR1) / RR1
    corr[RR_zero] = np.nan
    return corr, corr1


def angcor(nmin, nmax, flux):
    """
    :param nmin:
    :param nmax:
    :return:
    """
    fname_data = "/dataspace/sandeep/Angcor/TGSS_data/data_Cat/TGSS_%dmJy.txt" % flux
    ra_data = np.genfromtxt(fname_data,usecols=0, delimiter='\t')
    dec_data = np.genfromtxt(fname_data,usecols=1, delimiter='\t')

    data = np.asarray(ra_dec_to_xyz(ra_data, dec_data), order='F').T
    bins1 = 10 ** np.linspace(np.log10(0.001), np.log10(10), 50)
    bins_transform = angular_dist_to_euclidean_dist(bins1)

    for i in xrange(nmin, nmax):
        name = '/dataspace/sandeep/Angcor/TGSS_data/random_Cat/ranCat_CUTE_%d.txt' % i
        ra_R = np.genfromtxt(name, usecols=0, delimiter='\t')
        dec_R = np.genfromtxt(name, usecols=1, delimiter='\t')
        data_R = np.asarray(ra_dec_to_xyz(ra_R, dec_R), order='F').T
        xi, xi1 = two_point(data, bins_transform, method='standard', data_R=data_R, random_state=None)

        name = "/dataspace/sandeep/Angcor/TGSS_data/Corr_data/angCor_%dmJy_%d.txt" % \
               (flux, i)
        np.savetxt(name, zip(bins1, xi1, xi), delimiter=',', fmt='%0.6e', header='theta,Corr,AvgCorr')


if __name__ == "__main__":

    print "Enter which flux cut you want"
    cut = float(raw_input())

    """
    Cell_Count1 = Process(target=angcor, args=(1, 31, cut))
    Cell_Count1.start()
    Cell_Count2 = Process(target=angcor, args=(31, 61, cut))
    Cell_Count2.start()
    Cell_Count3 = Process(target=angcor, args=(61, 91, cut))
    Cell_Count3.start()
    Cell_Count4 = Process(target=angcor, args=(91, 121, cut))
    Cell_Count4.start()
    Cell_Count5 = Process(target=angcor, args=(121, 151, cut))
    Cell_Count5.start()
    Cell_Count6 = Process(target=angcor, args=(151, 181, cut))
    Cell_Count6.start()
    Cell_Count7 = Process(target=angcor, args=(181, 211, cut))
    Cell_Count7.start()
    Cell_Count8 = Process(target=angcor, args=(211, 241, cut))
    Cell_Count8.start()
    Cell_Count9 = Process(target=angcor, args=(241, 271, cut))
    Cell_Count9.start()
    Cell_Count10 = Process(target=angcor, args=(271, 301, cut))
    Cell_Count10.start()
    Cell_Count11 = Process(target=angcor, args=(301, 331, cut))
    Cell_Count11.start()
    Cell_Count12 = Process(target=angcor, args=(331, 361, cut))
    Cell_Count12.start()
    Cell_Count13 = Process(target=angcor, args=(361, 391, cut))
    Cell_Count13.start()
    Cell_Count14 = Process(target=angcor, args=(391, 421, cut))
    Cell_Count14.start()
    Cell_Count15 = Process(target=angcor, args=(421, 451, cut))
    Cell_Count15.start()
    Cell_Count16 = Process(target=angcor, args=(451, 481, cut))
    Cell_Count16.start()
    Cell_Count17 = Process(target=angcor, args=(511, 541, cut))
    Cell_Count17.start()
    Cell_Count18 = Process(target=angcor, args=(541, 571, cut))
    Cell_Count18.start()

    """
    Cell_Count19 = Process(target=angcor, args=(571, 601, cut))
    Cell_Count19.start()
    Cell_Count20 = Process(target=angcor, args=(601, 631, cut))
    Cell_Count20.start()
    Cell_Count21 = Process(target=angcor, args=(631, 661, cut))
    Cell_Count21.start()
    Cell_Count22 = Process(target=angcor, args=(661, 691, cut))
    Cell_Count22.start()
    Cell_Count23 = Process(target=angcor, args=(691, 721, cut))
    Cell_Count23.start()
    Cell_Count24 = Process(target=angcor, args=(721, 751, cut))
    Cell_Count24.start()
    Cell_Count25 = Process(target=angcor, args=(751, 781, cut))
    Cell_Count25.start()
    Cell_Count26 = Process(target=angcor, args=(781, 831, cut))
    Cell_Count26.start()
    Cell_Count27 = Process(target=angcor, args=(831, 861, cut))
    Cell_Count27.start()
    Cell_Count28 = Process(target=angcor, args=(861, 900, cut))
    Cell_Count28.start()

    """
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

    """
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



"""

 if sdensity:

        Dec_Bin = np.linspace(-53., 91,144)
        print Dec_Bin
        pixel_area =  hp.pixelfunc.nside2pixarea(NSIDE, degrees=True)

        dec=np.asarray(dec)
        pix_theta, pix_phi = 0.,0.

        SD = np.zeros(len(Dec_Bin), dtype=np.double)

        for i in xrange(len(Dec_Bin)):
            SUM=0
            if i< len(Dec_Bin)-1:
                index = (Dec_Bin[i]<= dec)*(dec<= Dec_Bin[i+1])
                for j in index:
                    if j == True:
                        SUM+=1
        SD[i]=SUM


        Dec_min = np.abs(np.subtract((-53.), 90.))
        Dec_Bin = np.arange(0., Dec_min+1.,1.)

        bin_area = np.zeros(len(Dec_Bin), dtype=np.double)
        print len(Dec_Bin)
        healpix_dec =np.zeros(npix, float)
        print npix

        for ipix in xrange(0, npix):

            pix_theta, pix_phi =hp.pixelfunc.pix2ang(NSIDE, ipix)
            healpix_dec[ipix] = pix_theta

        healpix_dec = np.degrees(healpix_dec)

        print len(bin_area)
        for i in xrange(len(Dec_Bin)):
            sum=0
            if i< len(Dec_Bin)-1:
                index = (Dec_Bin[i]<= healpix_dec)*(healpix_dec<= Dec_Bin[i+1])
                for j in index:
                    if j == True:
                        sum+=1
        bin_area[i]=sum


        bin_area*= pixel_area


        plt.figure(1, figsize=(8,6))
        plt.plot(Dec_Bin, SD/bin_area, 'b*')
        plt.yscale('log')
        plt.grid(True)
        plt.show()
"""









