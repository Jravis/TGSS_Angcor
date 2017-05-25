import numpy as np
from sklearn.neighbors import BallTree
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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

    corr = (factor ** 2 * DD - 2 * factor * DR + RR) / RR
    corr1 = (factor ** 2 * DD1 - 2 * factor * DR1 + RR1) / RR1
    #corr1 = (factor * DD/ DR) - 1.
    #corr2 = (factor**2*DD/ RR) - 1.

    corr[RR_zero] = np.nan

    return corr, corr1, DD, RR, DR 

"""
#****************************************************************
#Main

ra_32 = np.genfromtxt("dataCat.txt",usecols=0, delimiter='\t')
dec_32 =np.genfromtxt("dataCat.txt",usecols=1, delimiter='\t')


ra_100 = np.genfromtxt("dataCat_100.txt",usecols=0, delimiter='\t')
dec_100 =np.genfromtxt("dataCat_100.txt",usecols=1, delimiter='\t')


ra_300 = np.genfromtxt("dataCat_300.txt",usecols=0, delimiter='\t')
dec_300 =np.genfromtxt("dataCat_300.txt",usecols=1, delimiter='\t')



Nd1 = (1.0*len(ra_32))

data = np.asarray(ra_dec_to_xyz(ra_32, dec_32), order='F').T
#data1 = np.asarray(ra_dec_to_xyz(ra_100, dec_100), order='F').T
#data2 = np.asarray(ra_dec_to_xyz(ra_300, dec_300), order='F').T

bins1 = np.linspace(0.2,11,20)
bins2 = np.linspace(0.2,11,19)
bins_transform = angular_dist_to_euclidean_dist(bins1)
print bins1

bootstraps = []
bootstraps1 = []
bootstraps2 = []
seed = [ 55462255, 55555333, 56598313,57081391, 
        58343671, 58466453, 58947631, 59217501,
        59365582, 59534743]

Nbootstraps = 10

for i in xrange(Nbootstraps):

    name = "/home/sandeep/TGSS/Random_Catalog/Angcor/ranCat_%d.txt" % seed[i]
    
    ra_R = np.genfromtxt(name,usecols=0, delimiter='\t')
    dec_R =np.genfromtxt(name,usecols=1, delimiter='\t')
    Nr = (1.0*len(ra_R))
    data_R = np.asarray(ra_dec_to_xyz(ra_R, dec_R), order='F').T
    xi, xi1, dd, rr, dr = two_point(data, bins_transform, method='standard', data_R=data_R,
                               random_state=None)

#    xi1, dd1, rr1, dr1 = two_point(data1, bins_transform, method='standard', data_R=data_R,
#                               random_state=None)

#    xi2, dd2, rr2, dr2 = two_point(data2, bins_transform, method='standard', data_R=data_R,
#                               random_state=None)
    bootstraps.append(xi)
    bootstraps1.append(xi1)
#    bootstraps2.append(xi2)



bootstraps = np.asarray(bootstraps)
bootstraps1 = np.asarray(bootstraps1)
#bootstraps2 = np.asarray(bootstraps2)
print bootstraps.shape

Corr = np.mean(bootstraps, 0)
corr_err = np.std(bootstraps, 0, ddof=1)


Corr1 = np.mean(bootstraps1, 0)
corr_err1 = np.std(bootstraps1, 0, ddof=1)

#Corr2 = np.mean(bootstraps2, 0)
#corr_err2 = np.std(bootstraps2, 0, ddof=1)


#np.savetxt("bootstrapangCor.txt", zip(Corr, corr_err, Corr1, corr_err1, Corr2, corr_err2, bins2), delimiter='\t')
np.savetxt("bootstrap_BinnedSum.txt", zip(Corr, corr_err, Corr1, corr_err1, bins1, bins2), delimiter='\t')
#np.savetxt("RRDD.txt", zip(dd, rr, dr), delimiter='\t')

"""
Corr1 = np.genfromtxt("bootstrap_BinnedSum.txt", usecols=2)
corr_err1 = np.genfromtxt("bootstrap_BinnedSum.txt", usecols=3)

#Corr2 = np.genfromtxt("bootstrapangCor.txt", usecols=4)
#corr_err2 = np.genfromtxt("bootstrapangCor.txt", usecols=5)


Corr = np.genfromtxt("bootstrap_BinnedSum.txt", usecols=0)
corr_err = np.genfromtxt("bootstrap_BinnedSum.txt", usecols=1)


bins1 = np.linspace(0.2,11,20)
bins2 = np.linspace(0.2,11,19)

#w = []
#w1 = []
#for i in bins2:
#    w.append(0.5*(10**-1.5*i**(-1.6)+10**-1.63*i**(-1.2)))
x_l =[]
x_u =[]
y_l = []
y_u = []

for i in xrange(len(bins2)):
    if bins2[i] > 3.0:
        x_l.append(bins2[i])
        y_l.append(Corr[i])
    else:
        x_u.append(bins2[i])
        y_u.append(Corr[i])


print bins2
print max(x_u)
def fit_func(x1, a, b):
    return a*x1**-b

def fit_func1(x1, a, b, c):
    theta_0 = 3.0
    return a*x1**b + a*theta_0**(b-c)*x1**c


params = curve_fit(fit_func, x_u, y_u)

[a_u, b_u] = params[0]

params = curve_fit(fit_func, x_l, y_l)

[a_l, b_l] = params[0]


params = curve_fit(fit_func1, bins2, Corr)

[a, b, c] = params[0]


print a_l, a_u*3.0**(b_l-b_u)
x_u = np.asarray(x_u)
x_l = np.asarray(x_l)

w_u = fit_func(x_u, a_u, b_u)
w_l = fit_func(x_l, a_l, b_l)
w = fit_func1(bins2, a, b, c)

cos_u = 1-np.cos(np.radians(x_u)) 
sin_u = np.sin(np.radians(x_u)) 


cos_l = 1-np.cos(np.radians(x_l)) 
sin_l = np.sin(np.radians(x_l)) 

cos_bin = 1-np.cos(np.radians(bins2)) 
sin_bin = np.sin(np.radians(bins2)) 


dx_u = np.radians(x_u[1])-np.radians(x_u[0])
dx_l = np.radians(x_l[1])-np.radians(x_l[0])
y_u = np.multiply(cos_u, w_u)
y_l = np.multiply(cos_l, w_l)
dydx_u = np.gradient(y_u, dx_u)/sin_u
dydx_l = np.gradient(y_l, dx_l)/sin_l


y = np.multiply(cos_bin, Corr)
dx = np.radians(bins2[1])-np.radians(bins2[0])
dydx= np.gradient(y, dx)/sin_bin


plt.figure(1, figsize=(8,6))
plt.plot(bins2, Corr, 'b-', label='32mJy/beam',
        linewidth=3)
plt.plot(bins2, Corr1, linestyle='-', color='green', label='32mJy/beam',
        linewidth=2)

#plt.plot(bins2, w, linestyle='-', color='m',
#        linewidth=2)

plt.plot(x_u, w_u, '*-', color='m',
        linewidth=2)

plt.plot(x_l, w_l, '*-', color='m',
        linewidth=2)


plt.plot(x_u, dydx_u, '*-', color='r',
        linewidth=2)

plt.plot(x_l, dydx_l, '*-', color='k',
        linewidth=2)



plt.plot(bins2, dydx, '*-', color='y',
        linewidth=2)



plt.fill_between(bins2, Corr-corr_err, Corr+corr_err, alpha=0.5, edgecolor='gray', facecolor='gray')
plt.fill_between(bins2, Corr1-corr_err1, Corr1+corr_err1, alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
#plt.fill_between(bins2, Corr2-corr_err2, Corr2+corr_err2, alpha=0.5, edgecolor='y', facecolor='y')

plt.legend()
plt.grid(True)
plt.yscale("log")
plt.xscale("log")
plt.xlabel(r"$\theta$", fontsize=16)
plt.ylabel(r"W($\theta$)", fontsize=16)
plt.savefig("bootstrap_angCor.eps", dpi=100)


"""
plt.figure(5)
plt.plot(bins2, dd/(Nd1*(Nd1-1)*0.5),"bo-", label = 'DD')
plt.plot(bins2, rr/(Nr*(Nr-1)*0.5),"ro-", label='RR')
plt.plot(bins2, dr/(Nd1*Nr),"go-", label='DR')
plt.legend()
plt.grid(True)
plt.title("32mJy")
plt.savefig("angCor_Count1.eps", dpi=100)

plt.figure(6)
plt.plot(bins2, dd1/(Nd2*(Nd2-1)*0.5),"bo-", label = 'DD')
plt.plot(bins2, rr1/(Nr*(Nr-1)*0.5),"ro-", label='RR')
plt.plot(bins2, dr1/(Nd2*Nr),"go-", label='DR')
plt.legend()
plt.grid(True)
plt.title("100mJy")
plt.savefig("angCor_Count2.eps", dpi=100)


plt.figure(7)
plt.plot(bins2, dd2/(Nd3*(Nd3-1)*0.5),"bo-", label = 'DD')
plt.plot(bins2, rr2/(Nr*(Nr-1)*0.5),"ro-", label='RR')
plt.plot(bins2, dr2/(Nd3*Nr),"go-", label='DR')
plt.legend()
plt.grid(True)
plt.title("300mJy")
plt.savefig("angCor_Count3.eps", dpi=100)
"""
plt.show()
