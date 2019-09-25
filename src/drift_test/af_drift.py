#!/usr/bin/env python3

import numpy as np
from scipy.stats import binom, beta
from scipy.special import hyp2f1
import matplotlib.pyplot as plt


#mac = [0, 95]
#tac = [24, 172]
mac = [12, 168]
tac = [20, 172]

# Kimura (1964), Eq 4.10.
def kimura_drift(x, p, N, t, terms_for_convergence=50):
    phi = 0.0
    for i in range(1,terms_for_convergence):
        phi += p*(1-p)*i*(i+1)*(2*i+1)*hyp2f1(1-i, i+2, 2, p)*hyp2f1(1-i, i+2, 2, x)*np.exp(-i*(i+1)*t/(4*N))

    return phi

def af_pdf(k, n, x):
    """
    Alelle Frequency Probability Distribution Function.
    """
    return beta.pdf(x, 1+k, 1+(n-k))

def af_drift2(x, k, n, N, t):
    phi_p = np.zeros_like(x)
    dx = x[1] - x[0]
    for p, xx in zip(af_pdf(k, n, x), x):
        phi_p += kimura_drift(x, xx, N, t) * p
    return phi_p * dx

x = np.linspace(0, 1, num=1000)

plt.plot(x, af_pdf(mac[0], tac[0], x), "r", label="P1_MAF$\sim$Beta(1+{},1+{})".format(mac[0], tac[0]-mac[0]))
plt.plot(x, af_pdf(mac[1], tac[1], x), "b", label="P2_MAF$\sim$Beta(1+{},1+{})".format(mac[1], tac[1]-mac[1]))
plt.plot(x, af_pdf(sum(mac), sum(tac), x), "c", label="null (P1+P2) $\sim$Beta(1+{},1+{})".format(sum(mac), sum(tac)-sum(mac)))

af_mean = mac[0]/tac[0]
af_low, af_high = beta.interval(0.95, 1+mac[0], 1+tac[0]-mac[0])
p1_aflow_drift = kimura_drift(x, af_low, 1000, 50)
p1_afmean_drift = kimura_drift(x, af_mean, 1000, 50)
p1_afhigh_drift = kimura_drift(x, af_high, 1000, 50)
plt.plot(x, p1_aflow_drift, color="orange", label="drift(P1_MAF={:.3f}, N=1000, t=50)".format(af_low))
plt.plot(x, p1_afmean_drift, "m", label="drift(P1_MAF={:.3f}, N=1000, t=50)".format(af_mean))
plt.plot(x, p1_afhigh_drift, "k", label="drift(P1_MAF={:.3f}, N=1000, t=50)".format(af_high))

p1_af_drift = af_drift2(x, mac[0], tac[0], 1000, 50)
plt.plot(x, p1_af_drift, "g", lw=3, alpha=0.5, label="drift(P1_MAF$\sim$Beta(1+{},1+{}), N=1000, t=50)".format(mac[0], tac[0]-mac[0]))

plt.legend()
plt.xlabel("Minor Allele Frequency")
plt.ylabel("PDF")
#plt.title("")
plt.show()
