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
import collections
import numpy as np
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# SAM flag field
F_PAIRED        = 0x001 # the read is paired in sequencing
F_PAIR_MAPPED   = 0x002 # the read is mapped in a proper pair
F_UNMAPPED      = 0x004 # the query sequence itself is unmapped
F_MATE_UNMAPPED = 0x008 # the mate is unmapped
F_STRAND        = 0x010 # strand of the query (1 for reverse)
F_MATE_STRAND   = 0x020 # strand of the mate
F_FIRST_READ    = 0x040 # the read is the first read in a pair
F_SECOND_READ   = 0x080 # the read is the second read in a pair
F_SUBSEQUENT    = 0x100 # the alignment is not primary
F_QCFAIL        = 0x200 # QC failure
F_DUP           = 0x400 # optical or PCR duplicate
F_SUPP          = 0x800 # supplementary alignment

class ParseError(Exception):
    pass

def parse_sam(filename):
    lineno = 0
    with open(filename) as f:
        for line in f:
            lineno += 1
            if line[0] == '@':
                # ignore headers
                continue
            line = line.rstrip("\r\n")
            fields = line.split("\t")

            if len(fields) < 11:
                raise ParseError("{}: line {}: expected at least 11 columns, got {}.".format(filename, lineno, len(fields)))

            qname = fields[0]
            flag = int(fields[1])
            #rname = fields[2]
            #pos = fields[3]
            mapq = fields[4]
            #cigar = fields[5]
            #mrnm = fields[6]
            #mpos = fields[7]
            #isize = fields[8]
            seq = fields[9]
            #qual = fields[10]
            #opt = fields[11:]

            yield qname, flag, mapq, seq

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Plot fragment length histograms of un/mapped reads")
    parser.add_argument("--title", default=None, help="label for the title")
    parser.add_argument("--labels", help="comma separated labels for the input datasets")
    parser.add_argument("--max_len", type=int, default=10000, help="reads longer than this are ignored")
    parser.add_argument("in_sam", nargs='+', help="in1.sam [...inN.sam]")
    parser.add_argument("out_pdf", default="plot_fraglen_hist.pdf", help="out.pdf")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    lmin = 1e6
    lmax = 0

    labels = args.labels.split(",")
    if len(labels) != len(args.in_sam):
        print(file=sys.stderr)

    lengths = collections.defaultdict(lambda: np.zeros(args.max_len))

    for sam_file in args.in_sam:

        for _, flag, mapq, seq in parse_sam(sam_file):

            if flag & (F_PAIRED | F_UNMAPPED | F_SUBSEQUENT | F_QCFAIL | F_DUP | F_SUPP):
                    # Paired sequences are uninformative for length,
                    # ignore additional alignments and duplicates
                    continue

            length = len(seq)
            if length > args.max_len:
                continue

            lengths[sam_file][length] += 1

            if length < lmin:
                lmin = length
            elif length > lmax:
                lmax = length

    pdf = PdfPages(args.out_pdf)
    fig_w, fig_h = plt.figaspect(9.0/16.0)
    #fig_w, fig_h = plt.figaspect(3.0/4.0)

    fig1 = plt.figure(figsize=(fig_w, fig_h))
    ax1 = fig1.add_subplot(111)

    colours = ("#1b9e77", "#d95f02", "#7570b3")
#    colors = ("#66c2a5", "#fc8d62", "#8da0cb")

    bins =  np.arange(lmin, lmax+1) -0.5

    for sam_file, colour, label in zip(args.in_sam, colours, labels):
        f_lengths = lengths[sam_file][lmin:lmax+1] / sum(lengths[sam_file][lmin:lmax+1])
        ax1.bar(bins, f_lengths, width=0.5, alpha=0.5, color=colour, edgecolor=colour, label=label)

    ax1.set_xlabel("Length (bp)")
    ax1.set_ylabel("Frequency")
    ax1.legend()

    title = "Fragment length histogram"
    if args.title:
        title += " ({})".format(args.title)
    ax1.set_title(title)

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
