import numpy as np 
import warnings
from sklearn.neighbors import BallTree
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.pyplot import figure, show, rc


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
    """Convert ra & dec to Euclidean points

    Parameters
    ----------
    ra, dec : ndarrays

    Returns
    x, y, z : ndarrays
    """
    sin_ra = np.sin(ra * np.pi / 180.)
    cos_ra = np.cos(ra * np.pi / 180.)

    sin_dec = np.sin(np.pi / 2 - dec * np.pi / 180.)
    cos_dec = np.cos(np.pi / 2 - dec * np.pi / 180.)

    return (cos_ra * sin_dec,
            sin_ra * sin_dec,
            cos_dec)


def two_point(data, bins, method='standard', data_R=None, random_state=None):
    """Two-point correlation function

    Parameters
    ----------
    data : array_like
        input data, shape = [n_samples, n_features]
    bins : array_like
        bins within which to compute the 2-point correlation.
        shape = Nbins + 1
    method : string
        "standard" or "landy-szalay".
    data_R : array_like (optional)
        if specified, use this as the random comparison sample
    random_state : integer, np.random.RandomState, or None
        specify the random state to use for generating background

    Returns
    -------
    corr : ndarray
        the estimate of the correlation function within each bin
        shape = Nbins
    """
   # data = np.asarray(data)
   # bins = np.asarray(bins)

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
            counts_DD[i] = np.sum(BT_D.query_radius(data, bins[i],
                                                    count_only=True))
            counts_RR[i] = np.sum(BT_R.query_radius(data_R, bins[i],
                                                    count_only=True))

    DD1 = np.diff(counts_DD)
    RR1 = np.diff(counts_RR)

    DD = (counts_DD)
    RR = (counts_RR)


    # check for zero in the denominator
    RR_zero = (RR == 0)
    RR[RR_zero] = 1

   # if method == 'standard':
   #     corr = factor ** 2 * DD / RR - 1
   # elif method == 'landy-szalay':
    if sklearn_has_two_point:
        counts_DR = KDT_R.two_point_correlation(data, bins)
    else:
        counts_DR = np.zeros(Nbins + 1)
        for i in range(Nbins + 1):
            counts_DR[i] = np.sum(BT_R.query_radius(data, bins[i],
                                                        count_only=True))
    DR1 = np.diff(counts_DR)
    DR = (counts_DR)
    print "factor", factor
    corr = (factor ** 2 * DD - 2 * factor * DR + RR) / RR
    corr1 = (factor ** 2 * DD1 - 2 * factor * DR1 + RR1) / RR1
    #corr1 = (factor * DD/ DR) - 1.
    #corr2 = (DD*RR/ DR**2) - 1.

    corr[RR_zero] = np.nan

    return corr, corr1#, corr2, DD, RR, DR 

#****************************************************************
#Main

ra_32 = np.genfromtxt("../out_data/dataCat_4milijyPerBeam_8sigmaDete_Tflux.100Jy.txt",usecols=0, delimiter='\t')
dec_32 =np.genfromtxt("../out_data/dataCat_4milijyPerBeam_8sigmaDete_Tflux.100Jy.txt",usecols=1, delimiter='\t')


Nd1 = (1.0*len(ra_32))
data = np.asarray(ra_dec_to_xyz(ra_32, dec_32), order='F').T

bins1 = np.linspace(0.2,11,20)
bins2 = np.linspace(0.2,11,19)
bins_transform = angular_dist_to_euclidean_dist(bins1)
print bins1

seed = [ 55462255, 55555333, 56598313,57081391, 
        58343671, 58466453, 58947631, 59217501,
        59365582, 59534743]

Nbootstraps = 10
name = "../out_data/ranCat_CUTE_%d.txt" % seed[0]
print name
    
ra_R = np.genfromtxt(name, usecols=0, delimiter='\t')
dec_R =np.genfromtxt(name, usecols=1, delimiter='\t')
Nr = (1.0*len(ra_R))
data_R = np.asarray(ra_dec_to_xyz(ra_R, dec_R), order='F').T

xi, xi1 = two_point(data, bins_transform, method='standard', data_R=data_R,
               random_state=None)


bins2 = np.linspace(0.2, 11, 19)

costerm = 1-np.cos(np.radians(bins1)) 
sinterm = np.sin(np.radians(bins1)) 

dx = np.radians(bins1[1])-np.radians(bins1[0])
y = np.multiply(costerm, xi)
dydx = np.gradient(y, dx)/sinterm

print len(xi), len(xi1), len(dydx)
np.savetxt("angCor.txt", zip(xi, xi1, dydx, bins2), delimiter='\t')

plt.figure(1, figsize=(8,6))
plt.plot(bins2, xi1, "m*--", label='LS_Bindiff')
plt.plot(bins1, xi, "g*-", label='LS_Binsum')
plt.plot(bins1, dydx, "rD-", label='diffrential')
plt.xlabel(r"$\theta$")
plt.ylabel(r"$w(\theta)$")
plt.legend()
plt.grid(True)
plt.yscale("log")
plt.xscale("log")
#plt.savefig("angCor1.eps", dpi=100)

plt.show()



