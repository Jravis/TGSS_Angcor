import numpy as np
import matplotlib.pyplot as plt
import lmfit


name = "/dataspace/sandeep/Angcor/TGSS_data/Corr_data/50mJy_data/angCor_%dmJy_%d.txt" % (50, 1)
theta = np.genfromtxt(name, usecols=0, delimiter=',')
flx = [50, 60, 100, 200]
mean1 = np.zeros((4, len(theta)), dtype=np.float64)
std_dev = np.zeros((4, len(theta)), dtype=np.float64)

for flux in xrange(len(flx)):
    esti_bis_1 = np.zeros((899, 49), dtype=np.float64)

    for fn in xrange(1, 900):

        if flx[flux] == 50:
            if 481 <= fn < 511:
                name = "/dataspace/sandeep/Angcor/TGSS_data/Corr_data/%dmJy_data/angCor_%d.txt" % (flx[flux], fn)
            else:
                name = "/dataspace/sandeep/Angcor/TGSS_data/Corr_data/%dmJy_data/angCor_%dmJy_%d.txt" % (flx[flux], flx[flux], fn)
        else:
            name = "/dataspace/sandeep/Angcor/TGSS_data/Corr_data/%dmJy_data/angCor_%d.txt" % (flx[flux], fn)
        tmp = np.genfromtxt(name, usecols=1, delimiter=',')
        esti_bis_1[fn-1, :] = tmp

    mean1[flux, :] = np.mean(esti_bis_1, 0, dtype=np.float64)
    std_dev[flux, :] = np.std(esti_bis_1, 0, dtype=np.float64)


X = []
Y = []

y = mean1[0, :]
x = theta
y_err = std_dev[0, :]


index = (x >= 0.1) * (x <= 1.0)
data = np.zeros((3, len(x[index])), dtype=np.double)
data[0, :] = x[index]
data[1, :] = y[index]
data[2, :] = y_err[index]


pars = lmfit.Parameters()
pars.add_many(('a', 0.008), ('b', 0.5))


def residual(p):
    a = p['a'].value
    b = p['b'].value
    return (data[1, :] - (a*data[0, :]**(-b)))/data[2, :]

mini = lmfit.Minimizer(residual, pars)
#result = mini.minimize()
#print(lmfit.fit_report(result.params))

# first solve with Nelder-Mead
out1 = mini.minimize(method='Nelder')

# then solve with Levenberg-Marquardt using the
# Nelder-Mead solution as a starting point
out2 = mini.minimize(method='leastsq', params=out1.params)

lmfit.report_fit(out2.params, min_correl=0.5)

ci, trace = lmfit.conf_interval(mini, out2, sigmas=[0.68, 0.95],
                                trace=True, verbose=False)
print "FOR 0.1 <= theta <= 1.0) confidence interval and best fit"
print("--------------------------")
lmfit.printfuncs.report_ci(ci)

cx, cy, grid = lmfit.conf_interval2d(mini, out2, 'a', 'b', 30, 30)

#  ***********************************************************************

index = (x > 1.0) 
data = np.zeros((3,len(x[index])), dtype=np.double)
data[0, :] = x[index]
data[1, :] = y[index]
data[2, :] = y_err[index]


pars = lmfit.Parameters()
pars.add_many(('a', 0.009), ('b', 0.8))


def residual(p):
    a = p['a'].value
    b = p['b'].value
    return (data[1, :] - (a*data[0, :]**(-b)))/data[2, :]
mini = lmfit.Minimizer(residual, pars)

# first solve with Nelder-Mead
out1 = mini.minimize(method='Nelder')

# then solve with Levenberg-Marquardt using the
# Nelder-Mead solution as a starting point

out2 = mini.minimize(method='leastsq', params=out1.params)
lmfit.report_fit(out2.params, min_correl=0.5)
ci, trace = lmfit.conf_interval(mini, out2, sigmas=[0.68, 0.95],
                                trace=True, verbose=False)
print ""
print "FOR theta > 1.0) confidence interval and best fit"
print("--------------------------")
lmfit.printfuncs.report_ci(ci)
cx1, cy1, grid1 = lmfit.conf_interval2d(mini, out2, 'a', 'b', 30, 30)


plt.figure(1, figsize=(8, 6))
plt.contourf(cx, cy, grid, np.linspace(0.0, 1.0, 11), cmap=plt.cm.RdBu)
plt.xlabel('A', fontsize=16)
plt.colorbar()
plt.ylabel(r'$\gamma$', fontsize=16)
#plt.grid(True)
#plt.savefig("ParamContour.0.1-1deg.eps", dpi=100)

plt.savefig('/home/sandeep/TGSS_Fit_Plot/ParamContour.0.1-1deg.eps', dpi=100)

plt.figure(2, figsize=(8, 6))
plt.contourf(cx1, cy1, grid1, np.linspace(0.0,1.0,11), cmap=plt.cm.RdBu)
plt.xlabel('A', fontsize=16)
plt.colorbar()
plt.ylabel(r'$\gamma$', fontsize=16)
#plt.grid(True)
#plt.savefig("ParamContour.1degAbove.eps", dpi=100)

plt.savefig('/home/sandeep/TGSS_Fit_Plot/ParamContour.1degAbove.eps', dpi=100)

plt.figure(3)
plt.errorbar(x, y, yerr=y_err, fmt='.', color='crimson', ms=8, label='50mJy')
plt.plot(x, 0.00887917*x**-0.36569419, "g-", linewidth=2)

plt.fill_between(x, (0.00887917-0.00046)*x**-(0.36569419-0.06093),  (0.00887917+0.00047)*x**-(0.36569419+0.05908), alpha = 0.5, edgecolor='c',
                 facecolor='#AAAAAA')

plt.plot(x, 0.00891831*x**-0.78290750, "r-", linewidth=2)

plt.fill_between(x, (0.00891831-0.00046)*x**-(0.78290750-0.04496),  (0.00891831+0.00046)*x**-(0.78290750+0.04532), alpha = 0.5, edgecolor='c',
                 facecolor='skyblue')

plt.xscale('log')
plt.yscale('log')
plt.ylim(0.0006, 0.049)
plt.legend()
plt.xlim(0.09, 15)
plt.grid(which='both')
plt.legend()
plt.xlabel(r'$\theta$', fontsize='x-large', fontstyle='italic', weight='extra bold')
plt.ylabel(r'$\omega$', fontsize='x-large', fontstyle='italic', weight='extra bold')
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)

plt.savefig('/home/sandeep/TGSS_Fit_Plot/chi2_TGSS_fit.eps', dpi=100)
plt.show()




