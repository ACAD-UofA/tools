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
Produce chromosome coverage line plots for SNP data from VCFs.

e.g.
$ ./plot_chrom_cov_from_vcf.py \
        --bed /localscratch/BisonProjects/SNPs/9908_Bovid_SNP_10k.bed \
        --fai /localscratch/Refs/Bos_taurus/Bos_taurus_NoUnplaced_BWA6_2_2013_04/REFBTaurusNoUnplaced.fasta.fai \
        /localscratch/grg/150623_IMS_BLl_BHC04Bison10K/*BisonX*.vcf.gz
"""

from __future__ import print_function
import sys
import collections
from operator import itemgetter
import parsers

import numpy as np
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
#from matplotlib.ticker import FixedLocator, FixedFormatter
import matplotlib.ticker as ticker
from matplotlib.artist import setp

def chr_cmp(a, b):
    if a == b:
        return 0
    try:
        if int(a) < int(b):
            return -1
        else:
            return 1
    except ValueError:
        return cmp(a, b)

def parse_fai(filename):
    chrom_dict = {}
    with open(filename) as f:
        for line in f:
            line = line.rstrip("\r\n")
            fields = line.split("\t")

            chrom = fields[0]
            size = int(fields[1])

            if chrom.startswith("chr"):
                chrom = chrom[3:]
            
            chrom_dict[chrom] = size

    # sort by chrom
    return collections.OrderedDict(sorted(chrom_dict.items(), cmp=chr_cmp, key=itemgetter(0)))

def parse_bed(filename):
    # FIXME: merge into parsers.parse_bed()

    loci = collections.defaultdict(list)

    with open(filename) as f:
        for line in f:
            line = line.rstrip("\r\n")
            fields = line.split("\t")

            chrom = fields[0]
            #start = int(fields[1])
            end = int(fields[2])
            snpid = int(fields[3])

            if chrom.startswith("chr"):
                chrom = chrom[3:]

            loci[chrom].append(end)

    # sort by chrom
    loci = collections.OrderedDict(sorted(loci.items(), cmp=chr_cmp, key=itemgetter(0)))

    for chrom in loci:
        # sort by pos
        loci[chrom].sort()

    return loci

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Produce chromosome coverage line plots for SNP data from VCFs.")
    parser.add_argument("--bed", required=True, help="Bed intervals for SNPs.")
    parser.add_argument("--fai", required=True, help="Reference fasta index.")
    parser.add_argument("vcf", nargs="+", help="One or more VCFs (optionally gzipped).")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    chrom_dict = parse_fai(args.fai)
    loci = parse_bed(args.bed)

    total_dp = collections.defaultdict(lambda:collections.Counter())

    pdf = PdfPages("plot_chrom_cov_from_vcf.pdf")
    fig_w, fig_h = plt.figaspect(4.0/1.0)
    #fig_w, fig_h = plt.figaspect(9.0/16.0)
    #fig_w, fig_h = plt.figaspect(3.0/4.0)

    for filename in args.vcf:
        vparser = parsers.vcf_parser(filename)
        samples = next(vparser)

        x_pos = collections.defaultdict(lambda:collections.defaultdict(list))
        y_dp = collections.defaultdict(lambda:collections.defaultdict(list))

        max_dp = collections.defaultdict(int)

        for row in vparser:
            (lineno, line, chrom, pos, id, ref, alt, qual, filter, info_str, fmt_fields, gts) = row

            for sample in samples:
                dp = gts[sample].get("DP")
                if not dp:
                    print("Missing FORMAT/DP for chr{}:{}, sample ``{}'' at line {}".format(chrom, pos, sample, lineno), file=sys.stderr)
                    sys.exit(1)

                dp = int(dp)

                if dp > max_dp[sample]:
                    max_dp[sample] = dp

                x_pos[sample][chrom].append(pos)
                y_dp[sample][chrom].append(dp)

                total_dp[chrom][pos] += dp

        max_chr_len = max(chrom_dict.values())

        for sample in samples:
            fig = plt.figure(figsize=(fig_w*5, fig_h*5))
            fig.suptitle(sample, fontsize=32)

            cols=1
            gs = gridspec.GridSpec(len(chrom_dict)/cols + len(chrom_dict)%cols, cols)

            ax = None
            for (chrom, chrom_len), ss in zip(chrom_dict.iteritems(), gs):
                ax = fig.add_subplot(ss, sharex=ax)

                ax.set_title("chr{}".format(chrom))
                ax.xaxis.set_tick_params(direction="out")
                #ax.xaxis.set_major_formatter(ticker.NullFormatter()) #lambda x, _: "{}Mb".format(int(x/1e6))))
                setp(ax.get_xticklabels(), visible=False)
                ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: str(int(x))))

                # FIXME: set axis bounds to extent of chromosome
                (x, y, w, h) = ax.get_position().bounds
                new_w = w * float(chrom_len)/max_chr_len
                ax.set_position([x, y, new_w, h])

                try:
                    ax.vlines(loci[chrom], 0, 0.2, 'r')#, lw=0.3)
                except KeyError:
                    pass

                ax.vlines(x_pos[sample][chrom], 0, y_dp[sample][chrom], 'k')#, lw=0.3)
                ax.set_ylim(top=max_dp[sample])
                #ax.yaxis.set_major_locator(ticker.LinearLocator(max_dp[sample]))

                ax.vlines(chrom_len, 0, max_dp[sample], 'b', lw=1)

            ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: "{}Mb".format(int(x/1e6))))
            setp(ax.get_xticklabels(), visible=True)
            ax.set_xlim(1, max_chr_len)

            gs.tight_layout(fig, rect=[0, 0, 1, 0.97])
            pdf.savefig(figure=fig)


    total_max_dp = 0
    for chrom in chrom_dict:
        try:
            tmp_dp = max(total_dp[chrom].values())
        except ValueError:
            continue
        if tmp_dp > total_max_dp:
            total_max_dp = tmp_dp


    fig = plt.figure(figsize=(fig_w*5, fig_h*5))
    fig.suptitle("Total", fontsize=32)

    cols=1
    gs = gridspec.GridSpec(len(chrom_dict)/cols + len(chrom_dict)%cols, cols)

    ax = None
    for (chrom, chrom_len), ss in zip(chrom_dict.iteritems(), gs):
        ax = fig.add_subplot(ss, sharex=ax)

        ax.set_title("chr{}".format(chrom))
        ax.xaxis.set_tick_params(direction="out")
        #ax.xaxis.set_major_formatter(ticker.NullFormatter()) #lambda x, _: "{}Mb".format(int(x/1e6))))
        setp(ax.get_xticklabels(), visible=False)
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: str(int(x))))

        # FIXME: set axis bounds to extent of chromosome
        #(x, y, w, h) = ax.get_position().bounds
        #new_w = w * float(chrom_len)/max_chr_len
        #ax.set_position([x, y, new_w, h])

        xx = [a[1] for a in sorted(total_dp[chrom].items(), key=itemgetter(0))]

        try:
            ax.vlines(loci[chrom], 0, 0.2, 'r')#, lw=0.3)
            ax.vlines(loci[chrom], 0, xx, 'k')#, lw=0.3)
        except KeyError:
            pass

        ax.set_ylim(top=total_max_dp)
        #ax.yaxis.set_major_locator(ticker.LinearLocator(total_max_dp))

        ax.vlines(chrom_len, 0, total_max_dp, 'b', lw=1)

    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: "{}Mb".format(int(x/1e6))))
    setp(ax.get_xticklabels(), visible=True)
    ax.set_xlim(1, max_chr_len)

    gs.tight_layout(fig, rect=[0, 0, 1, 0.97])
    pdf.savefig(figure=fig)

    pdf.close()
