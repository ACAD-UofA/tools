#!/usr/bin/env python

from __future__ import print_function
import sys
import collections
import os.path
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
from matplotlib import cm
import numpy as np
from scipy import stats

def parse_methcounts(filename, dpmax, countmin):
    mc = collections.defaultdict(list)
    with open(filename) as f:
        next(f)
        for lineno, line in enumerate(f,2):
            try:
                line = line.rstrip()
                if not line:
                    continue
                m, u, n = map(int, line.split("\t"))
            except:
                print("{}:``{}''".format(lineno, line), file=sys.stderr)
                raise
            if m+u > dpmax or n < countmin:
                continue
            mc[m+u].append((m,u,n))
    return mc

def methdist_f(x, methcounts, f):
    yy = np.empty((len(methcounts),len(x)))
    nsites = 0
    for i, (m,u,n) in enumerate(methcounts):
        yy[i,:] = n*f(x, m+1, u+1)
        nsites += n
    return np.sum(yy, axis=0) / nsites

methdist_pdf = lambda a,b: methdist_f(a,b,stats.beta.pdf)
methdist_cdf = lambda a,b: methdist_f(a,b,stats.beta.cdf)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: {} methcounts.txt methdist.pdf".format(sys.argv[0]),
                file=sys.stderr)
        exit(1)

    dpmin = 2
    dpmax = 15 #1024
    countmin = 0

    title = os.path.basename(sys.argv[1])
    if title.endswith(".txt"):
        title = title[:-4]

    mc = parse_methcounts(sys.argv[1], dpmax, countmin)
#    mc_all = [x for mcx in mc.values() for x in mcx if x[0]+x[1]<dpmax]

    mc_all = []
    methprops = []
    depth = []
    for mcx in mc.values():
        mc_all.extend(mcx)
        for m,u,n in mcx:
            methprops.extend(n*[float(m)/(m+u)])
            depth.extend(n*[m+u])

    pdf = PdfPages(sys.argv[2])
    fig_w, fig_h = plt.figaspect(3.0/4.0)
    fig1 = plt.figure(figsize=(fig_w, 3*fig_h))
    #fig1 = plt.figure(figsize=(8.27,11.69)) # a4 portrait
    #fig1 = plt.figure(figsize=(11.69,8.27)) # a4 landscape
    gs1 = gridspec.GridSpec(3, 1)
    ax1 = fig1.add_subplot(gs1[0])
    ax2 = fig1.add_subplot(gs1[1])
    ax3 = fig1.add_subplot(gs1[2], projection='3d')

    ax1.set_title(title)

    ax1.hist(depth, normed=True)
    ax1.set_xlabel("depth")
    ax1.set_ylabel("frequency")

    npoints = 100
    x0 = np.linspace(0, 1-1.0/npoints, npoints)
    x1 = x0+1.0/npoints
    x = x0+0.5/npoints
    y = np.arange(dpmin,dpmax)
    z = methdist_pdf(x, mc_all)

    ax2.plot(x, z, color="black", lw=2)
    ax2.hist(methprops, bins=10, color="red", normed=True)
    ax2.set_xlabel("methylation proportion")
    ax2.set_ylabel("frequency")

    X,Y = np.meshgrid(x, y)
    Z = np.empty((len(y),len(x)))
    for dp in y:
        Z[dp-dpmin,:] = methdist_pdf(x, mc[dp])

    surf = ax3.plot_surface(X,Y,Z, cmap=cm.coolwarm, rstride=5, cstride=20)
    ax3.set_xlabel("methylation proportion")
    ax3.set_ylabel("depth")
    ax3.set_zlabel("frequency")

    fig1.colorbar(surf, shrink=0.5, aspect=5)

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
