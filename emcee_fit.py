import lmfit
import corner

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

"""
index = (x >= 0.1) * (x <= 1.0)
data = np.zeros((3, len(x[index])), dtype=np.double)
data[0, :] = x[index]
data[1, :] = y[index]
data[2, :] = y_err[index]


pars = lmfit.Parameters()
pars.add_many(('A', 0.008), ('b', 0.5))


def residual(p1):
    a = p1['A'].value
    b = p1['b'].value
    return data[1, :] - (a*data[0, :]**(-b))  #+ c*data[0, :]**(-d))


def lnprob(p):
    resid = residual(p)
    s = data[2, :]
    resid *= 1 / s
    resid *= resid
    resid += np.log(2 * np.pi * s ** 2)
    return -0.5 * np.sum(resid)

mi = lmfit.Minimizer(residual, pars)
# first solve with Nelder-Mead
out1 = mi.minimize(method='Nelder')


print("Statistics for 0.1<=theta<=1.0")
mini = lmfit.Minimizer(lnprob, out1.params)
res = mini.emcee(burn=300, steps=1000, thin=10, params=mi.params)

print("median of posterior probability distribution")
print('------------------------------------------')
lmfit.report_fit(res.params)

fig1 = corner.corner(res.flatchain, labels=res.var_names, truths=list(res.params.valuesdict().values()),
                    truth_color='skyblue')

fig1.savefig('/home/sandeep/emcee_TGSS_fit_1.eps', dpi=100)
highest_prob = np.argmax(res.lnprob)
hp_loc = np.unravel_index(highest_prob, res.lnprob.shape)
mle_soln = res.chain[hp_loc]
for i, par in enumerate(pars):
    pars[par].value = mle_soln[i]

print("\nMaximum likelihood Estimation")
print('-----------------------------')
print(pars)
"""

print""
print("Statistics for theta>1.0")
index = (x > 1.0)  #* (x <= 1.0)
data = np.zeros((3, len(x[index])), dtype=np.double)
data[0, :] = x[index]
data[1, :] = y[index]
data[2, :] = y_err[index]


pars = lmfit.Parameters()
pars.add('A', value=0.008, min=0.001, max=0.2)
pars.add('b', value=0.75)#, min=0.7, max=0.8)


def residual(p1):
    a = p1['A'].value
    b = p1['b'].value
    return data[1, :] - (a*data[0, :]**(-b))


def lnprob(p):
    resid = residual(p)
    s = data[2, :]
    resid *= 1 / s
    resid *= resid
    resid += np.log(2 * np.pi * s ** 2)
    return -0.5 * np.sum(resid)

mi = lmfit.Minimizer(residual, pars)
# first solve with Nelder-Mead
out1 = mi.minimize(method='Nelder')
mini = lmfit.Minimizer(lnprob, out1.params)

res = mini.emcee(burn=400, steps=6000, thin=10, nwalkers=100,  is_weighted=False, params=mi.params, seed=1230127)

print('------------------------------------------')
print("median of posterior probability distribution")
print('------------------------------------------')

lmfit.report_fit(res.params)


fig1 = corner.corner(res.flatchain, bins=300, labels=res.var_names,
                    truths=list(res.params.valuesdict().values()),title_fmt=".4f",
                    quantiles=[0.16, 0.5, 0.84], show_titles=True, labels_args={"fontsize": 40})

fig1.set_size_inches(8, 6)
fig1.savefig('/dataspace/sandeep/Angcor/TGSS_data/Corr_data/emcee_TGSS_fit1.png', dpi=600)


highest_prob = np.argmax(res.lnprob)
hp_loc = np.unravel_index(highest_prob, res.lnprob.shape)
mle_soln = res.chain[hp_loc]
for i, par in enumerate(pars):
    pars[par].value = mle_soln[i]

print("\nMaximum likelihood Estimation")
print('-----------------------------')
print(pars)

plt.show()