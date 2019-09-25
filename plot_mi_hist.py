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

from __future__ import print_function
import sys
import gzip
import numpy as np
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class ParseError(Exception):
    pass

def parse_tsv2dictlist(filename, keep_cols):

    if filename.endswith(".gz"):
        xopen = gzip.open
    else:
        xopen = open

    dlist = []
    lineno = 1

    with xopen(filename) as f:
        hline = next(f)
        headers = hline.rstrip().split("\t")
        for line in f:
            lineno += 1
            line = line.rstrip()
            fields = line.split("\t")

            if len(headers) != len(fields):
                raise ParseError("{}: line {}: number of fields ({}) does not match the number of headers ({}).".format(filename, lineno, len(fields), len(headers)))

            d = {}
            for key, val in zip(headers, fields):
                if key in keep_cols:
                    d[key] = val
            dlist.append(d)

    return headers, dlist

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: {} beast_out.log[.gz]".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    tiso_cols = ("tmrca(NewWorld_A)", "tmrca(NewWorld_B)", "tmrca(NewWorld_C)", "tmrca(NewWorld_D1)", "tmrca(NewWorld_D4h3a)")
    tdiv_cols = ("tmrca(A)", "tmrca(B)", "tmrca(C)", "tmrca(D)")

    keep_cols = set()
    keep_cols.update(tiso_cols)
    keep_cols.update(tdiv_cols)

    headers, dlist = parse_tsv2dictlist(sys.argv[1], keep_cols)

    # time of complete isolation, time of divergence
    tiso = np.empty(len(dlist))
    tdiv = np.empty(len(dlist))
    tiso_sum = 0
    tdiv_sum = 0

    for i, d in enumerate(dlist):
        tiso[i] = min(float(d[key]) for key in tiso_cols)
        tdiv[i] = max(float(d[key]) for key in tdiv_cols)
        tiso_sum += tiso[i]
        tdiv_sum += tdiv[i]

    tiso_ci_low, tiso_median, tiso_ci_high = np.percentile(tiso, [2.5, 50, 97.5])
    tdiv_ci_low, tdiv_median, tdiv_ci_high = np.percentile(tdiv, [2.5, 50, 97.5])
    tiso_mean = float(tiso_sum) / len(dlist)
    tdiv_mean = float(tdiv_sum) / len(dlist)

    print("\tTime of complete isolation\tTime of divergence")
    print("mean\t{}\t{}".format(tiso_mean, tdiv_mean))
    print("median\t{}\t{}".format(tiso_median, tdiv_median))
    print("95% CI lower bound\t{}\t{}".format(tiso_ci_low, tdiv_ci_low))
    print("95% CI upper bound\t{}\t{}".format(tiso_ci_high, tdiv_ci_high))

    pdf = PdfPages("plot_mi_hist.pdf")
    fig_w, fig_h = plt.figaspect(9.0/16.0)
    fig1 = plt.figure(figsize=(fig_w, fig_h))
    ax1 = fig1.add_subplot(111)

    n, bins, patches = ax1.hist(tiso, bins=200, normed=True)
    ax1.set_xlabel("Time of complete isolation")
    ax1.set_ylabel("Frequency")
    ax1.set_title("Histogram of Isolation Times")
    for b, p in zip(bins, patches):
        if b < tiso_ci_low or b > tiso_ci_high:
            p.set_color('r')

    ax1.axvline(12472.5281921177, 0, 1, color='black')
    ax1.axvline(18362.3383578888, 0, 1, color='black')

    plt.tight_layout()
    pdf.savefig()

    fig1 = plt.figure(figsize=(fig_w, fig_h))
    ax1 = fig1.add_subplot(111)

    n, bins, patches = ax1.hist(tdiv, bins=200, normed=True)
    ax1.set_xlabel("Time of divergence")
    ax1.set_ylabel("Frequency")
    ax1.set_title("Histogram of Divergence Times")
    for b, p in zip(bins, patches):
        if b < tdiv_ci_low or b > tdiv_ci_high:
            p.set_color('r')
    ax1.axvline(21778.9945493312, 0, 1, color='black')
    ax1.axvline(38714.9132320423, 0, 1, color='black')

    plt.tight_layout()
    pdf.savefig()

    pdf.close()
