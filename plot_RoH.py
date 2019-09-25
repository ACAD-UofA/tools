#!/usr/bin/env python

from __future__ import print_function
import sys
import numpy as np
from collections import defaultdict

import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

class SortError(Exception):
    pass

def parse_geno(fn):
    with open(fn) as f:
        for line in f:
            yield line.strip()

def parse_ind(fn):
    with open(fn) as f:
        for line in f:
            fields = line.split()
            sample = fields[0]
            group = fields[2]
            yield sample, group

def parse_snp(fn):
    with open(fn) as f:
        for line in f:
            fields = line.split()
            #id = fields[0]
            chrom = fields[1]
            #gpos = fields[2]
            pos = int(fields[3])
            #ref = fields[4]
            #alt = fields[5]
            yield chrom, pos

def agg_fixedwin(generator, winsz):
    last_chrom = None
    data_list = []
    i = 0
    for entry in generator:
        chrom, pos, data = entry[0], entry[1], entry[2:]
        if last_chrom != chrom or pos > i+winsz:
            if last_chrom != chrom:
                i = 0
                last_chrom = chrom
            while pos > i+winsz:
                yield last_chrom, i, i+winsz,  data_list
                i += winsz
                data_list = []
        elif pos < i:
            raise SortError("{}:{} follows {}:{}".format(chrom, pos, chrom, last_pos))

        data_list.append(data)

    yield last_chrom, i, i+winsz, data_list

def flatter(xxx):
    for xx in xxx:
        for x in xx:
            yield x

flatten = lambda g: list(flatter(g))

def agg_sliding(agg, nchunks):
    data_q = []
    pos_q = []
    last_chrom = None

    for chrom, pos0, pos1, data_list in agg:
        if last_chrom != chrom:
            last_chrom = chrom
            data_q = []
            pos_q = []

        data_q.append(data_list)
        pos_q.append((pos0, pos1))

        if len(data_q) < nchunks:
            continue

        yield chrom, pos_q[0][0], pos_q[-1][1], flatten(data_q)

        data_q.pop(0)
        pos_q.pop(0)

def parse_args():
    import argparse

    def km_int(s):
        if s[-1] in "kK":
            s2 = s[:-1]
            m = 1000
        elif s[-1] in "mM":
            s2 = s[:-1]
            m = 1000000
        else:
            s2 = s
            m = 1
        try:
            x = int(s2)
        except:
            raise argparse.ArgumentTypeError("unrecognised size `{}'".format(s))
        return x*m

    parser = argparse.ArgumentParser(description="plot Runs of Homozygosity (RoH) from SNP data")
    parser.add_argument("-w", "--window", type=km_int, default="5m", help="use windows of size WINDOW [%(default)s]")
    parser.add_argument("-s", "--step", type=km_int, default="200k", help="slide windows with step size STEP [%(default)s]")
    parser.add_argument("-m", "--missing-threshold", type=float, default=1.0/3.0, help="exclude windows with more missing data than this proportion [%(default)s]")
    parser.add_argument("--scale", type=float, default=1.0, help="scale the plot [%(default)s]")
    parser.add_argument("--aspect", type=float, default=9.0/16.0, help="aspect ratio of plots [%(default)s]")
    parser.add_argument("data_pfx", type=str, help="path to the EIGENSOFT format .geno .snp and .ind files")
    parser.add_argument("out_pdf", type=str, help="path for output pdf")
    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = parse_args()

    geno_fn = args.data_pfx + ".geno"
    snp_fn = args.data_pfx + ".snp"
    ind_fn = args.data_pfx + ".ind"

    indlist = list(parse_ind(ind_fn))
    n_inds = len(indlist)

    Hsum = np.zeros(n_inds, dtype=int)
    Nsum = np.zeros(n_inds, dtype=int)
    Msum = np.zeros(n_inds, dtype=int)

    chrom_list = []
    chrom_set = set()

    last_chrom = None
    chrom_bounds = []
    cum_pos = 0
    last_pos = 0
    roh2d = defaultdict(list)
    X_all = [1]
    Xc = defaultdict(lambda:[1])

    get_gts = ((c,p,g) for (c,p),g in zip(parse_snp(snp_fn), parse_geno(geno_fn)))
    nchunks = args.window / args.step
    aggregator = agg_sliding(agg_fixedwin(get_gts, args.step), nchunks)

    for chrom, pos0, pos1, gts_list in aggregator:
        midpos = int(round(float(pos0 + pos1 + 1)/2.0))

        if chrom not in chrom_set:
            chrom_list.append(chrom)
            chrom_set.add(chrom)
        if chrom != last_chrom and last_chrom is not None:
            chrom_bounds.append(cum_pos)
            last_pos = cum_pos
        cum_pos = last_pos + midpos
        X_all.append(cum_pos)
        last_chrom = chrom

        H = np.zeros(n_inds, dtype=int)
        N = np.zeros(n_inds, dtype=int)
        M = np.zeros(n_inds, dtype=int)

        for gts in gts_list:
            for i, gt in enumerate(gts[0]):
                if gt == "9":
                    M[i] += 1
                else:
                    N[i] += 1
                    if gt == "1":
                        H[i] += 1

        roh = np.empty(n_inds, dtype=float)
        for i, (h,n,m) in enumerate(zip(H,N,M)):
            if n < 10 or float(m)/(m+n) > args.missing_threshold:
                roh[i] = float('nan')
            else:
                roh[i] = float(n-h)/n
            Hsum[i] += h
            Nsum[i] += n
            Msum[i] += m

        roh2d[chrom].append(roh)
        Xc[chrom].append(midpos)


    roh2d_all = []
    for c in chrom_list:
        roh2d_all.extend(roh2d[c])
    roh2d["Autosome"] = roh2d_all
    Xc["Autosome"] = X_all
    chrom_list = ["Autosome"] + chrom_list

    print("sample", "group", "HET%", "MISS%", sep="\t")
    for (ind,group), H, N, M in zip(indlist, Hsum, Nsum, Msum):
        print(ind, group, "{:.2f}".format(100.0*H/N),
                          "{:.2f}".format(100.0*M/(M+N)), sep="\t")

    chrom_bounds2 = [0]+chrom_bounds+[cum_pos]
    chrom_mid = []
    for i in range(1,len(chrom_bounds2)):
        mid = chrom_bounds2[i-1] + (chrom_bounds2[i] - chrom_bounds2[i-1]) / 2.0
        chrom_mid.append(mid)

    pdf = PdfPages(args.out_pdf)

    Y = np.arange(n_inds+1)
    #cmap = matplotlib.cm.get_cmap('coolwarm')
    cmap = matplotlib.cm.get_cmap('bwr')
    cmap.set_bad(color='black')

    for chrom in chrom_list:

        fig_w, fig_h = plt.figaspect(args.aspect)
        fig1 = plt.figure(figsize=(args.scale*fig_w, args.scale*fig_h))
        gs1 = gridspec.GridSpec(1, 1)
        ax1 = fig1.add_subplot(gs1[0])

        X = Xc[chrom]
        Z = np.ma.array(roh2d[chrom], mask=np.isnan(roh2d[chrom])).T
        heatmap = ax1.pcolormesh(X, Y, Z, cmap=cmap, edgecolors='None', vmin=0.9, vmax=Z.max(), antialiased=False, rasterized=True)
        heatmap._is_stroked = False

        if chrom == "Autosome":
            ax1.set_xticks(chrom_mid)
            ax1.set_xticklabels(np.arange(1,len(chrom_mid)+1), size="6")
            ax1.tick_params(axis='x', which='major', length=0)

            ax1.tick_params(axis='x', which='minor', width=1, length=4)
            ax1.set_xticks(chrom_bounds, minor=True)
        else:
            xticks = np.arange(2e7, 3e8, 2e7)
            xlabels = [str(int(x/1e6))+" mbp" for x in xticks]
            ax1.set_xticks(xticks)
            ax1.set_xticklabels(xlabels, size="6")

        ax1.set_title("Homozygosity ({}{})".format("Chr" if chrom!="Autosome" else "", chrom))

        ax1.set_yticks(Y[:-1]+0.5)
        ytl = ax1.set_yticklabels([ind for ind,_ in indlist], size="4")

        ax1.set_xlim(X[0]-0.5, X[-1]+0.5)
        ax1.set_ylim(Y[0], Y[-1])


        plt.colorbar(heatmap, orientation="vertical", ax=ax1)

        plt.tight_layout()
        pdf.savefig()
        plt.close()
    pdf.close()
