#!/usr/bin/env python

import sys
import os.path
import collections
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def parse_lengths(filename):
    hist = []
    fits = collections.defaultdict(dict)
    with open(filename) as f:
        for line in f:
            line = line.rstrip("\r\n")
            if len(line) == 0:
                continue
            if line[0] == "#":
                fields = line[1:].split(";")
                for ff in fields[1:]:
                    k, v = ff.split("=")
                    fits[fields[0]][k] = float(v)
                continue
            hist.append(int(line))
    return np.array(hist, dtype=float), fits

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Plot fragment length histogram")
    parser.add_argument("-m", "--maxlen", type=int, default=0, help="maximum length to plot")
    parser.add_argument("--lognorm", action="store_true", default=False, help="plot lognorm")
    parser.add_argument("--exp2", action="store_true", default=False, help="plot lognorm")
    parser.add_argument("--exp2mix", action="store_true", default=False, help="plot lognorm")
    parser.add_argument("--lognormmix", action="store_true", default=False, help="plot lognormmix")
    parser.add_argument("--weibull", action="store_true", default=False, help="plot weibull")
    parser.add_argument("--rayleigh", action="store_true", default=False, help="plot rayleigh")
    parser.add_argument("-t", "--title", help="title")
    parser.add_argument("in_hist", help="in.hist")
    parser.add_argument("out_pdf", default="plot_length_distrib.pdf", help="out.pdf")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    hist, fits = parse_lengths(args.in_hist)

    if args.maxlen:
        hist = hist[:args.maxlen]

    pdf = PdfPages(args.out_pdf)
    #fig_w, fig_h = plt.figaspect(9.0/16.0)
    fig_w, fig_h = plt.figaspect(3.0/4.0)
    fig1 = plt.figure(figsize=(fig_w, fig_h))
    ax1 = fig1.add_subplot(111)

    area = np.sum(hist)
    #hist /= area

    # black, orange, sky blue, blueish green, yellow, dark blue, vermillion, reddish purple
#    pallete = ("#000000", "#906000", "#357090", "#006050", "#959025", "#004570", "#804000", "#806070")
    # blue/green, orange, purple
    pallete = ("#000000", "#1b9e77", "#d95f02", "#7570b3")
    bins =  np.arange(1, len(hist)+1)
    width=0.8
    ax1.bar(bins-width/2.0, hist, width=width, color=pallete[0], edgecolor=pallete[0])

    show_legend = False

    if args.lognorm:
        meanlog = fits["lnorm"]["meanlog"]
        sdlog = fits["lnorm"]["sdlog"]
        fit_ln = stats.lognorm.pdf(bins, sdlog, scale=np.exp(meanlog))
        ax1.plot(bins, fit_ln*area, color=pallete[3], ls="--", lw=3, label="lognormal")
        show_legend = True

    if args.exp2:
        a = fits["exp2"]["a"]
        b = fits["exp2"]["b"]
        fit_exp2 = np.exp(-a*bins-b/bins)
        ax1.plot(bins, fit_exp2*area, color=pallete[2], ls="-", lw=3, label="exp2")
        show_legend = True

    if args.exp2mix:
        a = fits["exp2mix"]["a"]
        b = fits["exp2mix"]["b"]
        c = fits["exp2mix"]["c"]
        y = fits["exp2mix"]["y"]
        ssqe = fits["exp2mix"]["SSqE"]
        fit_e1 = y*exp(-a*bins -b/bins)
        fit_e2 = (1-y)*exp(-c*bins -b/bins)
        fit_emix = fit_e1 + fit_e2
        ax1.plot(bins, fit_emix*area, color=pallete[2], ls="-", lw=3, label="lognormal1+lognormal2")
        ax1.plot(bins, fit_e1*area, color=pallete[1], ls=":", lw=3, label="lognormal1")
        ax1.plot(bins, fit_e2*area, color=pallete[1], ls=":", lw=3, label="lognormal2")
        show_legend = True

    if args.lognormmix:
        meanlog1 = fits["lnormmix"]["meanlog1"]
        sdlog1 = fits["lnormmix"]["sdlog1"]
        meanlog2 = fits["lnormmix"]["meanlog2"]
        sdlog2 = fits["lnormmix"]["sdlog2"]
        y = fits["lnormmix"]["gamma"]
        ssqe = fits["lnormmix"]["SSqE"]
        fit_ln1 = y*stats.lognorm.pdf(bins, sdlog1, scale=np.exp(meanlog1))
        fit_ln2 = (1-y)*stats.lognorm.pdf(bins, sdlog2, scale=np.exp(meanlog2))
        fit_ln = fit_ln1 + fit_ln2
        ax1.plot(bins, fit_ln*area, color=pallete[2], ls="-", lw=3, label="lognormal1+lognormal2")
        ax1.plot(bins, fit_ln1*area, color=pallete[1], ls=":", lw=3, label="lognormal1")
        ax1.plot(bins, fit_ln2*area, color=pallete[1], ls=":", lw=3, label="lognormal2")
        show_legend = True

#   meanlog = fits["lnormgeom"]["meanlog"]
#   sdlog = fits["lnormgeom"]["sdlog"]
#   p = fits["lnormgeom"]["p"]
#   #print(meanlog, sdlog, p)
#    fit_lng = stats.lognorm.pdf(bins, sdlog, scale=np.exp(meanlog)) * stats.geom.pmf(bins, p)
#   fit_lng = np.zeros(len(bins))
#   for i in range(len(fit_lng)):
#       ff = 0
#       for j in range(i,i+100):
#           if j == len(fit_ln):
#               break
#           #ff += fit_ln[j] * (1-p)**(j-i) * p
#           ff += fit_ln[j] * stats.poisson.pmf(j-i, p)
#       fit_lng[i] = ff
#    ax1.plot(bins, fit_lng*area, color=pallete[2], ls="--", lw=3, label="lnormgeom")

    if args.weibull:
        delta = fits["weibull"]["delta"]
        eta = fits["weibull"]["eta"]
        d0 = fits["weibull"]["d0"]
        ssqe = fits["weibull"]["SSqE"]
        fit_wb = delta/eta * ((bins-d0)/eta)**(delta-1) * np.exp(-((bins-d0)/eta)**delta)
        ax1.plot(bins, fit_wb*area, color=pallete[2], ls=":", lw=3, label="weibull")
        show_legend = True

    if args.rayleigh:
        loc = fits["rayleigh"]["loc"]
        scale = fits["rayleigh"]["scale"]
        ssqe = fits["rayleigh"]["SSqE"]
        fit_rl = stats.rayleigh.pdf(bins, loc=loc, scale=scale)
        ax1.plot(bins, fit_rl*area, color=pallete[1], ls="-.", lw=3, label="rayleigh")
        show_legend = True

    ax1.xaxis.set_tick_params(direction="inout")
    ax1.set_xlabel("Length (bp)")
    ax1.set_ylabel("Frequency")
    if show_legend:
        ax1.legend(numpoints=1)

    if args.title:
        ax1.set_title(args.title)

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
