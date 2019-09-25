#!/usr/bin/env python

from __future__ import print_function
import sys
import collections
import operator
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.stats import linregress

def parse_isize(fn, contigs, min_reads=500, skip=1):
    l_counts = collections.defaultdict(lambda:collections.defaultdict(dict))
    gc_counts = collections.defaultdict(lambda:collections.defaultdict(dict))
    with open(fn) as f:
        while skip:
            next(f)
            skip -= 1
        for line in f:
            fields = line.split()
            sample = fields[0]
            chrom = fields[1]
            length = int(fields[2])
            nreads = int(fields[3])
            gc = int(fields[4])
            l_counts[sample][chrom][length] = nreads
            gc_counts[sample][chrom][length] = gc

    # pull out samples with more than min_reads reads for every chromosome
    samples = []
    for s in l_counts.keys():
        n_list = [sum(l_counts[s][ch].values()) for ch in contigs]
        if min(n_list) > min_reads:
            samples.append(s)

    mean_lists_by_s = collections.defaultdict(list)
    sd_lists_by_s = collections.defaultdict(list)
    mean_lists_by_ch = collections.defaultdict(list)
    sd_lists_by_ch = collections.defaultdict(list)

    for ch in contigs:
        for s in samples:
            n = sum(l_counts[s][ch].values())

            # mean read length
            m_sum = 0
            for l, nr in l_counts[s][ch].iteritems():
                m_sum += l*nr
            mean = float(m_sum) / n

            # variance of read length
            sum_sq = 0
            for l, nr in l_counts[s][ch].iteritems():
                sum_sq += (l-mean)*(l-mean)*nr
            var = float(sum_sq) / (n-1)
            sd = np.sqrt(var)

            mean_lists_by_s[s].append(mean)
            mean_lists_by_ch[ch].append(mean)
            sd_lists_by_s[s].append(sd)
            sd_lists_by_ch[ch].append(sd)

    return (mean_lists_by_ch, sd_lists_by_ch,
            np.array([mean_lists_by_s[s] for s in samples]),
            np.array([sd_lists_by_s[s] for s in samples]),
            samples)

def parse_fai(fn):
    chrom_dict = {}
    with open(fn) as f:
        for line in f:
            fields = line.split()
            chrom = fields[0]
            size = int(fields[1])
            if chrom.startswith("chr"):
                chrom = chrom[3:]
            chrom_dict[chrom] = size
    return chrom_dict

def plot_regress(ax, x, y, title=None, colour="black"):

    ax.scatter(x, y, facecolor=colour, edgecolors=colour, lw=0.5, s=60)

    slope, intercept, r_val, p_val, std_err = linregress(x, y)
    min_x, max_x = np.min(x), np.max(x)
    min_x = min_x - 0.01*abs(min_x)
    max_x = max_x + 0.01*abs(max_x)

    lx = np.array([min_x, max_x])
    ly = slope*lx + intercept
    ax.plot(lx, ly, '--', color=colour)

    if title is None:
        title = ""
    else:
        title += ": "
    ax.set_title('{}n={}, $r^2$={:.3f}, $p$={:.3f}'.format(title, len(x), r_val*r_val, p_val))

    #ylim = ax.get_ylim()
    #ax.set_ylim([min(0, np.min(y)), ylim[1]])
    ax.set_xlim([min_x, max_x])

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: {} isize.txt fasta.fai".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    isize_fn = sys.argv[1]
    fai_fn = sys.argv[2]

    nchr = 29
    contigs = map(str, range(1,nchr+1)) # + ["X"]
    mean, sd, mean_rows, sd_rows, samples = parse_isize(isize_fn, contigs)

    chr_dict = parse_fai(fai_fn)
    chr_sizes = np.array([chr_dict[ch] for ch in contigs])

    pdf = PdfPages("isize.pdf")
    aspect = 9.0/16.0
    #aspect = 3.0/4.0
    fig_w, fig_h = plt.figaspect(aspect)
    fig1 = plt.figure(figsize=(2*fig_w, 2*fig_h))
    gs1 = gridspec.GridSpec(4, 1)
    ax1 = fig1.add_subplot(gs1[0])
    ax2 = fig1.add_subplot(gs1[1], sharex=ax1)
    ax3 = fig1.add_subplot(gs1[2])
    ax4 = fig1.add_subplot(gs1[3])
    #ax3 = fig1.add_subplot(gs1[2:4], sharex=ax1)
    #ax4 = fig1.add_subplot(gs1[4:6], sharex=ax1)

    x = np.arange(len(contigs))
    colour = "black"

    ax1.violinplot([np.log(mean[ch]) for ch in contigs], x, points=60, widths=0.7, showmedians=True, bw_method=0.5)
    ax1.set_ylabel("log(Mean read length)")
    ax1.set_xticks(x)
    ax1.set_xticklabels(contigs)
    ax1.set_xlim([min(x)-0.5, max(x)+0.5])

    hbar_mean = np.log(np.median([np.median(mean[ch]) for ch in contigs]))
    ax1.hlines(hbar_mean, min(x)-0.5, max(x)+0.5, linestyle=':', color=colour)

    ax2.violinplot([np.log(sd[ch]) for ch in contigs], x, points=60, widths=0.7, showmedians=True, bw_method=0.5)
    ax2.set_ylabel("log(StDev read length)")

    hbar_sd = np.log(np.median([np.median(sd[ch]) for ch in contigs]))
    ax2.hlines(hbar_sd, min(x)-0.5, max(x)+0.5, linestyle=':', color=colour)

    #ax3.bar(x, [len(mean[ch]) for ch in contigs])
    #ax4.bar(x, [len(sd[ch]) for ch in contigs])


    pvals = []
    rvals = []
    for row in mean_rows:
        slope, intercept, r_val, p_val, std_err = linregress(chr_sizes, row)
        pvals.append(p_val)
        if p_val < 0.05:
            rvals.append(r_val)

    ax3.hist(pvals, bins=20)
    ax3.set_title("histogram of $p$-vals (regression of mean read length against chromosome size)")
    ax3.set_xlabel("$p$-val")
    ax4.hist(rvals, bins=20)
    ax4.set_title("histogram of $r$-vals (regression of mean read length against chromosome size)")
    ax4.set_xlabel("$r$-val (Pearson correlation)")

    plt.tight_layout()
    pdf.savefig()
    plt.close()


    with pdf:
        for s, m_r, sd_r in zip(samples, mean_rows, sd_rows):
            fig_w, fig_h = plt.figaspect(aspect)
            fig1 = plt.figure(figsize=(2*fig_w, 2*fig_h))
            gs1 = gridspec.GridSpec(2, 1)
            ax1 = fig1.add_subplot(gs1[0])
            ax2 = fig1.add_subplot(gs1[1], sharex=ax1)

            plot_regress(ax1, np.log(chr_sizes), np.log(m_r), "{}: chrsize vs. mean read length".format(s))
            ax1.set_xlabel("log(chromosome length)")
            ax1.set_ylabel("log(mean read length)")
            plot_regress(ax2, np.log(chr_sizes), np.log(sd_r), "{}: chrsize vs. SD read length".format(s))
            ax2.set_xlabel("log(chromosome length)")
            ax2.set_ylabel("log(SD of read length)")

            plt.tight_layout()
            pdf.savefig()
            plt.close()

