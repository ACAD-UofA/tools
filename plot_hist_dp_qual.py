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
Plot histograms of one or more vcf's QUAL and sample DP.

e.g.
$ ./plot_hist_dp_qual.py /localscratch/jsoubrier/10k/BisonX/run_raxml/*.vcf.gz
"""

from __future__ import print_function
import sys
import collections
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import parsers

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("usage: {} in1.vcf [in2.vcf...inN.vcf]".format(sys.argv[0]))
        sys.exit(1)

    pdf = PdfPages("plot_hist_dp_qual.pdf")
    fig_w, fig_h = plt.figaspect(9.0/16.0)
    #fig_w, fig_h = plt.figaspect(3.0/4.0)

    for vcf_in in sys.argv[1:]:

        vparser = parsers.vcf_parser(vcf_in)
        samples = next(vparser)
        if (len(samples) != 1):
            print("Error: {}: multisample vcf is unexpected.".format(vcf_in), file=sys.stderr)
            sys.exit(1)

        sample = samples[0]

        qual_list = []
        dp_list = []
        gt_count = collections.Counter()

        for row in vparser:
            (lineno, line, chrom, pos, id, ref, alt, qual, filter, info_str, fmt_fields, gts) = row

            if qual == ".":
                print("Error: {}: line {}: missing QUAL field.".format(vcf_in, lineno), file=sys.stderr)
                sys.exit(1)

            qual_list.append(float(qual))

            for field in ("DP", "GT"):
                if field not in gts[sample]:
                    print("Error: {}: line {}: missing {} field for sample {}.".format(vcf_in, lineno, field, sample), file=sys.stderr)
                    sys.exit(1)

            dp_list.append(int(gts[sample]["DP"]))

            gt_str = gts[sample].get("GT")
            if gt_str == "./.":
                gt = 0
            else:
                gt = sum([int(x) for x in gt_str.split("/")])
            gt_count[gt] += 1

        print("gt_count", gt_count)

        fig1 = plt.figure(figsize=(fig_w, fig_h))

        gs1 = gridspec.GridSpec(1, 2)
        ax1, ax2 = [fig1.add_subplot(ss) for ss in gs1]

        ax1.hist(qual_list, bins=100)
        del qual_list
        ax1.set_xlabel("QUAL")
        ax1.set_ylabel("Frequency")

        ax2.hist(dp_list, bins=30)
        del dp_list
        ax2.set_xlabel("DP")

        fig1.suptitle(sample, fontsize=10)

        fig1.text(0.5, 0.95, "Total SNPs: {}, REF/REF: {}, REF/ALT: {}, ALT/ALT: {}".format(sum(gt_count.values()), gt_count[0], gt_count[1], gt_count[2]), ha="center", va="top", fontsize=10)
        del gt_count

        #gs1.tight_layout(fig1)
        pdf.savefig()#bbox_extra_artists=(leg1,))

    pdf.close()
