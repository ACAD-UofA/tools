#!/usr/bin/env python
# Copyright (c) 2015 Graham Gower <graham.gower@gmail.com>
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

"""
Plot histograms of how many samples cover a site.

e.g.
$ ./plot_samples_covered.py \
        --bed /localscratch/BisonProjects/SNPs/9908_Bovid_SNP_10k.bed \
        /localscratch/grg/150623_IMS_BLl_BHC04Bison10K/*BisonX*.vcf.gz
"""

from __future__ import print_function
import sys
import collections
import parsers

import numpy as np
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Plot histograms of how many samples cover SNP sites, from VCFs.")
    parser.add_argument("--bed", required=True, help="Bed intervals for SNPs.")
    parser.add_argument("vcf", nargs="+", help="One or more VCFs (optionally gzipped).")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    samples_covered = collections.Counter()
    _, loci_list = parsers.parse_bed(args.bed)
    for chrom, pos, _ in loci_list:
        samples_covered[(chrom, pos)] = 0
    all_samples = []

    pdf = PdfPages("plot_samples_covered.pdf")
    #fig_w, fig_h = plt.figaspect(4.0/1.0)
    fig_w, fig_h = plt.figaspect(9.0/16.0)
    #fig_w, fig_h = plt.figaspect(3.0/4.0)

    for filename in args.vcf:
        vparser = parsers.vcf_parser(filename)
        samples = next(vparser)
        all_samples.extend(samples)

        for row in vparser:
            (lineno, line, chrom, pos, id, ref, alt, qual, filter, info_str, fmt_fields, gts) = row

            for sample in samples:
                #dp_str = gts[sample].get("DP")
                #if not dp_str:
                #    print("Error: {}: line {}: Missing FORMAT/DP for sample ``{}'' at chr{}:{}.".format(filename, lineno, sample, chrom, pos), file=sys.stderr)
                #    sys.exit(1)
                #dp = int(dp_str)

                samples_covered[(chrom, pos)] += 1

    hist_bins = range(0, len(all_samples)+2)

    fig = plt.figure(figsize=(fig_w, fig_h))
    #fig.suptitle("") #, fontsize=10)

    gs = gridspec.GridSpec(2, 1)
    ax1 = None
    ax1, ax2 = [fig.add_subplot(ss, sharex=ax1) for ss in gs]

    hh, _, _ = ax1.hist(samples_covered.values(), bins=hist_bins, align="left")
    ax1.set_xlabel("$N$ samples with $DP\geq{}1$")
    ax1.set_ylabel("Number of loci")
    ax1.set_ylim([0, 1.2*max(hh)])
    for x, h in zip(hist_bins, hh):
        ax1.text(x, h+0.05*max(hh), str(int(h)), color="grey", ha="center", weight="bold")

    hh, _, _ = ax2.hist(samples_covered.values(), bins=hist_bins, align="left", cumulative=-1)
    ax2.set_xlabel("$N$ or more samples with $DP\geq{}1$")
    ax2.set_ylabel("Number of loci")
    ax2.set_ylim([0, 1.2*max(hh)])
    for x, h in zip(hist_bins, hh):
        ax2.text(x, h+0.05*max(hh), str(int(h)), color="grey", ha="center", weight="bold")

    for ax in (ax1, ax2):
        ax.set_xticks(hist_bins[:-1])
        ax.set_xlim([-0.5, len(all_samples)+0.5])

    gs.tight_layout(fig)#, rect=[0, 0, 1, 0.95])
    pdf.savefig(figure=fig)

    pdf.close()
