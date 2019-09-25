#!/usr/bin/env python

#GEN_TIME=8.0
#Ne=1e4
#TM0_LOWER=(1e5/GEN_TIME)/(2.0*Ne)
#TM1_UPPER=(5.9e6/GEN_TIME)/(2.0*Ne)
#TDIFF=(TM1_UPPER-TM0_LOWER)

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

pdf = PdfPages("chain.pdf")
fig_w, fig_h = plt.figaspect(9.0/16.0)
#fig_w, fig_h = plt.figaspect(3.0/4.0)
fig1 = plt.figure(figsize=(fig_w, fig_h))
gs1 = gridspec.GridSpec(2, 2)
ax1, ax2, ax3, ax4 = [fig1.add_subplot(ss) for ss in gs1]

next(sys.stdin) # skip header
y, _, t = zip(*[(float(f) for f in line.split("\t")[0:3]) for line in sys.stdin])
t = np.array(t)
y = np.array(y)
burnin=1000
chain_len=10000
t_chain = t[:chain_len]
y_chain = y[:chain_len]
t_burnt = t[burnin:]
y_burnt = y[burnin:]

# histogram weights, to normalise plots as PMFs
weights = np.ones_like(t_burnt)/float(len(t_burnt))

#ax1.plot(y_chain, linestyle="none", marker='o', ms=1)
ax1.plot(y_chain, 'k', ms=0)
#ax1.set_ylim((0.85, 1.0))
ax1.set_title("MCMC")
ax1.set_ylabel("$\gamma$")
ax1.set_xlabel("position in chain")

ax2.set_title("Stationary distribution ($\gamma$)")
ax2.hist(y_burnt, normed=False, weights=weights, bins=100)
ax2.set_xlabel("$\gamma$ (proportion of (AB,C) vs. (A,BC) gene trees)")

#ax3.plot(t_chain, linestyle="none", marker='o', ms=1)
ax3.plot(t_chain, 'k', ms=0)
#ax3.set_ylim((0.4, 0.8))
ax3.set_title("MCMC")
ax3.set_ylabel("$t$")
ax3.set_xlabel("position in chain")

ax4.set_title("Stationary distribution ($t$)")
ax4.hist(t_burnt, normed=False, weights=weights, bins=100)
ax4.set_xlabel("$t$ (time between hybridisation and TMRCA, $2N_e$ generations)")

ax1.set_xlim([0, len(y_chain)])
ax3.set_xlim([0, len(t_chain)])

plt.tight_layout()
pdf.savefig()
pdf.close()


pdf = PdfPages("y_vs_t.pdf")
#fig_w, fig_h = plt.figaspect(9.0/16.0)
fig_w, fig_h = plt.figaspect(3.0/4.0)
fig1 = plt.figure(figsize=(fig_w, fig_h))
ax1 = fig1.add_subplot(111)

_,_,_,im = ax1.hist2d(y_burnt, t_burnt, normed=False, weights=weights, bins=100)
#ax1.plot(y[1000:], t[1000:], linestyle="none", marker='o', ms=1)
ax1.set_xlabel("$\gamma$ (proportion of (AB,C) vs. (A,BC) gene trees)")
ax1.set_ylabel("$t$ ($2N_e$ generations)")
fig1.colorbar(im)

plt.tight_layout()
pdf.savefig()
pdf.close()
