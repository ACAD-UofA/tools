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
Scree plot from eigensoft .eval file.
"""

from __future__ import print_function
import sys
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import numpy as np

def parse_pca_evec(filename):
    evec = []
    with open(filename) as f:
        for line in f:
            evec.append(float(line.strip()))
    return evec

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: {} eigout.eval plot_out.pdf".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    evec = parse_pca_evec(sys.argv[1])

    pdf = PdfPages(sys.argv[2])
    fig_w, fig_h = plt.figaspect(9.0/16.0)
    #fig_w, fig_h = plt.figaspect(3.0/4.0)
    fig1 = plt.figure(figsize=(fig_w, fig_h))
    ax1 = fig1.add_subplot(111)

    components = np.arange(1, len(evec)+1)
    ax1.bar(components-0.4, evec)
    #ax1.set_xticks(np.arange(5, len(evec)+1, 5))
    ax1.set_xticks([])
    ax1.set_ylabel("Eigenvalue")
    ax1.set_xlabel("Component")

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
