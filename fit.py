import warnings
import numpy as np
from sklearn.neighbors import BallTree
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.pyplot import figure, show
from scipy.optimize import fmin
from scipy.optimize import curve_fit
import lmfit
#matplotlib.use('WXAgg')

def chi2(pars, data1):
    x, y, err = data1[0,:], data1[1,:], data1[2,:]
    residue = (y-model(pars, x))
    residue*= residue
    return residue.sum()


name = '/home/sandeep/CUTE-master/CUTE_angcor_data/NSIDE_1024_Corr_data/60miliJy/CUTE_Corr_55462255.dat'
x = np.genfromtxt(name, usecols=0)
y = np.genfromtxt(name, usecols=1)
e = np.genfromtxt(name, usecols=2)
dd = np.genfromtxt(name, usecols=3)

X = []
Y = []
for i in xrange(len(x)):
    if x[i] >= 0.2:
        X.append(x[i])
        Y.append(y[i])

np.savetxt("fit_test.txt",zip(X,Y), delimiter='\t')

index = (x >= 0.10) * (x<=1.0)
data = np.zeros((3,len(x[index])), dtype= np.double)
data[0,:]= x[index]
data[1,:]= y[index]
data[2,:]= e[index]


pars = lmfit.Parameters()
pars.add_many(('a', 0.008), ('b', 0.5))

def residual(p):
    a = p['a'].value
    b = p['b'].value
    return ( data[1,:] - (a*data[0,:]**(-b)) )

mini = lmfit.Minimizer(residual, pars)
#result = mini.minimize()
#print(lmfit.fit_report(result.params))

# first solve with Nelder-Mead
out1 = mini.minimize(method='Nelder')

# then solve with Levenberg-Marquardt using the
# Nelder-Mead solution as a starting point
out2 = mini.minimize(method='leastsq', params=out1.params)


lmfit.report_fit(out2.params, min_correl=0.5)

ci, trace = lmfit.conf_interval(mini, out2, sigmas=[0.68,0.95],
                                        trace=True, verbose=False)
lmfit.printfuncs.report_ci(ci)

cx, cy, grid = lmfit.conf_interval2d(mini, out2, 'a','b',30,30)
#***********************************************************************
index = (x > 1.0) 
data = np.zeros((3,len(x[index])), dtype= np.double)
data[0,:]= x[index]
data[1,:]= y[index]
data[2,:]= e[index]


pars = lmfit.Parameters()
pars.add_many(('a', 0.009), ('b', 0.8))

def residual(p):
    a = p['a'].value
    b = p['b'].value
    return ( data[1,:] - (a*data[0,:]**(-b)) )

mini = lmfit.Minimizer(residual, pars)

# first solve with Nelder-Mead
out1 = mini.minimize(method='Nelder')

# then solve with Levenberg-Marquardt using the
# Nelder-Mead solution as a starting point

out2 = mini.minimize(method='leastsq', params=out1.params)
lmfit.report_fit(out2.params, min_correl=0.5)
ci, trace = lmfit.conf_interval(mini, out2, sigmas=[0.68,0.95],
                                        trace=True, verbose=False)
lmfit.printfuncs.report_ci(ci)
cx1, cy1, grid1 = lmfit.conf_interval2d(mini, out2, 'a','b',30,30)


#*********************************************************************
plt.figure(1, figsize=(8,6))
#plt.contourf(cx, cy, grid, levels=[0.0,0.68,0.95], cmap=plt.cm.viridis)
plt.contourf(cx, cy, grid, np.linspace(0.0,1.0,11), cmap=plt.cm.viridis)
plt.xlabel('A', fontsize=16)
plt.colorbar()
plt.ylabel(r'$\gamma$', fontsize=16)
#plt.grid(True)
#plt.savefig("ParamContour.0.1-1deg.eps", dpi=100)

plt.figure(2, figsize=(8,6))
#plt.contourf(cx1, cy1, grid1, levels=[0.0,0.68,0.95], cmap=plt.cm.viridis)
plt.contourf(cx1, cy1, grid1, np.linspace(0.0,1.0,11), cmap=plt.cm.viridis)
plt.xlabel('A', fontsize=16)
plt.colorbar()
plt.ylabel(r'$\gamma$', fontsize=16)
#plt.grid(True)
#plt.savefig("ParamContour.1degAbove.eps", dpi=100)


#plt.figure(3)
#plt.plot(x,y,"b*", ms = 10)
#plt.plot(x,0.00989660*x**-0.33876366,"r-", linewidth=2)
#plt.plot(x,0.00995607*x**-0.732417,"g-", linewidth=2)
#plt.xscale('log')
#plt.yscale('log')
#plt.ylim(0.0006, 0.045)
#plt.xlim(0.05, 15 )
#plt.xlabel(r'$\theta$')
#plt.ylabel(r'\omega')
#plt.grid(True)

plt.figure(4, figsize=(8,6))
plt.plot(x,dd,"r*-", ms = 10)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\theta$', fontsize=16)
plt.ylabel('DD', fontsize=16)
plt.grid(True)
#plt.savefig("DataPairSepration.eps", dpi=100)

plt.show()




