#!/usr/bin/env python

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

pdf = PdfPages("ll.pdf")
#fig_w, fig_h = plt.figaspect(9.0/16.0)
fig_w, fig_h = plt.figaspect(3.0/4.0)
fig1 = plt.figure(figsize=(2*fig_w, 2*fig_h))
gs1 = gridspec.GridSpec(2, 2)
axes = [fig1.add_subplot(ss, projection='3d') for ss in gs1]

with open("llik.txt") as f:
    next(f) # skip header
    y, t, ab, bc, ac, data = zip(*((float(v) for v in line.rstrip("\r\n").split("\t")) for line in f))
labels = ("(AB,C)", "(A,BC)", "(B,AC)", "Data")
col4 = ("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")
col3 = ("#1b9e77", "#d95f02", "#7570b3")

for ax, label, ll, col in zip(axes, labels, (ab, bc, ac, data), col4):
    #ax.scatter(y, t, ll, c=col, edgecolor=col, s=1, label=label)
    ax.plot_trisurf(y, t, ll, cmap=cm.jet, linewidth=0.0)
    ax.set_title(label)
    ax.set_xlabel("$gamma$")
    ax.set_ylabel("$t$")
    ax.set_zlabel("$log(likelihood)$")

#ax1.set_xlim([-0.05, 1.05])
#ax2.set_xlim([-0.05, 1.05])

#plt.tight_layout()
pdf.savefig()
pdf.close()
