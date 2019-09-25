#!/usr/bin/env python

from scipy.stats import beta, poisson
import numpy as np

# age distribution, with lower limit 13yrs, mean 29yrs.
agedist = beta(2, 5, loc=13, scale=58)

mt_mu = 10**-5 # mito mutations/site/generation
mt_sz = 16500 # mito length

def count_mutations(c_time, niter=1000):

    mutations = np.zeros(niter)

    for i in xrange(niter):
        l1, l2 = 0, 0 # lineage ages
        ng1, ng2 = 0, 0 # number of generations per lineage
        while 1:
            age = agedist.rvs(1)
            if l1+age > c_time:
                break
            l1 += age
            ng1 += 1
        while 1:
            age = agedist.rvs(1)
            if l2+age > c_time:
                break
            l2 += age
            ng2 += 1

        # sum of poisson r.v.s are poisson with sum of means
        mutations[i] = poisson.rvs(1, mt_mu*mt_sz*(ng1+ng2))

    return mutations

times = np.arange(1e3, 1e4, 1e3)
mu_obs = np.zeros(len(times))

for i, tm in enumerate(times):
    mutations = count_mutations(tm)
    n_generations = float(tm)/agedist.mean()
    mu_obs[i] = np.mean(mutations) / (2.0*mt_sz*n_generations)
    print("tm={}".format(tm), mu_obs[i])

import matplotlib.pyplot as plt
plt.plot(times, mu_obs, 'kx')
plt.xlabel("Time of divergence (years)")
plt.ylabel("$\mu$ (mutations/site/generation)")
plt.show()

