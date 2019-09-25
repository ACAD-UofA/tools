#!/usr/bin/env python

from __future__ import print_function
import sys
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import numpy as np
import math
import operator

def setCustomHatchWidth(customWidth):
    """
    Monkeypatch the pdf writing for hatches, to change the hatch line width
    """

    # make sure you have the correct imports,
    # they may differ depending on the matplotlib version
    import matplotlib.backends.backend_pdf
    from matplotlib.externals import six
    from matplotlib.backends.backend_pdf import Name, Op
    from matplotlib.transforms import Affine2D
    from matplotlib.path import Path

    def _writeHatches(self):
        hatchDict = dict()
        sidelen = 72.0
        for hatch_style, name in six.iteritems(self.hatchPatterns):
            ob = self.reserveObject('hatch pattern')
            hatchDict[name] = ob
            res = {'Procsets':
                   [Name(x) for x in "PDF Text ImageB ImageC ImageI".split()]}
            self.beginStream(
                ob.id, None,
                {'Type': Name('Pattern'),
                 'PatternType': 1, 'PaintType': 1, 'TilingType': 1,
                 'BBox': [0, 0, sidelen, sidelen],
                 'XStep': sidelen, 'YStep': sidelen,
                 'Resources': res})

            stroke_rgb, fill_rgb, path = hatch_style
            self.output(stroke_rgb[0], stroke_rgb[1], stroke_rgb[2],
                        Op.setrgb_stroke)
            if fill_rgb is not None:
                self.output(fill_rgb[0], fill_rgb[1], fill_rgb[2],
                            Op.setrgb_nonstroke,
                            0, 0, sidelen, sidelen, Op.rectangle,
                            Op.fill)

            self.output(customWidth, Op.setlinewidth)

            # TODO: We could make this dpi-dependent, but that would be
            # an API change
            self.output(*self.pathOperations(
                Path.hatch(path),
                Affine2D().scale(sidelen),
                simplify=False))
            self.output(Op.stroke)

            self.endStream()
        self.writeObject(self.hatchObject, hatchDict)

    matplotlib.backends.backend_pdf.PdfFile.writeHatches = _writeHatches

def parse_hets(filename):
    h1 = []
    n1 = []
    with open(filename) as f:
        next(f) # skip header
        for line in f:
            line = line.rstrip()
            fields = line.split("\t")
            assert(len(fields) == 3)
            _, n, h = map(int, fields)
            if n == 0:
                break
            h1.append(float(h)/n)
            n1.append(n)

    return h1, n1

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: {} out.pdf hets.txt".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    setCustomHatchWidth(1)

    pdf = PdfPages(sys.argv[1])
    #fig_w, fig_h = plt.figaspect(9.0/16.0)
    fig_w, fig_h = plt.figaspect(3.0/4.0)
    fig1 = plt.figure(figsize=(fig_w, fig_h))

    gs1 = gridspec.GridSpec(2, 1)
    ax1 = fig1.add_subplot(gs1[0])
    ax2 = fig1.add_subplot(gs1[1], sharex=ax1)
    #ax1 = fig1.add_subplot(111)

# black, orange, sky blue, blueish green, yellow, dark blue, vermillion, reddish purple
#    pallete = ("#000000", "#906000", "#357090", "#006050", "#959025", "#004570", "#804000", "#806070")

    #pal = ["#1b9e77", "#d95f02", "#7570b3"]
    pal = ["#006050", "#806070", "#959025"]
    facecols = ['lightblue', '#660000', 'yellow']
    edgecols = ['darkblue', '#660000', 'darkyellow']
    vlinestyles = ['-', '--', '-x']
    fills = [True, False, True]
    hatches = ['', 'xx', '--']
    alpha = 1.0 #0.6
    barwidth = 1.0
    baroffset = 0.5
    xlim = 0
    ylim = False

    for col, ec, hatch, fill, vls, fn in zip(facecols, edgecols, hatches, fills, vlinestyles, sys.argv[2:]):
        h, n = parse_hets(fn)
        n100 = np.sum(n[1:])
        n05 = 0.05*n100
        n95 = 0.95*n100
        n99 = 0.99*n100
        sum_n = 0
        edf_n = [0.0,]
        i05 = False
        i95 = False
        i99 = False
        for i, nn in enumerate(n[1:], start=2):
            sum_n += nn
            edf_n.append(float(sum_n)/n100)
            if not i05 and sum_n > n05:
                i05 = i
            if not i95 and sum_n > n95:
                i95 = i
            if not i99 and sum_n > n99:
                i99 = i
        #print(i05, i95, i99)

        dph = list(enumerate(h[i05-1:i99], start=i05))
        dph = sorted(dph, key=operator.itemgetter(1))
        sum_n = 0
        for _dp, _h in dph:
            sum_n += n[_dp-1]
            edf = float(sum_n)/n100
            print(_dp, _h, edf, sep="\t")


        x = np.arange(1, len(h)+1)
        ax1.bar(x-baroffset, edf_n, width=barwidth, alpha=alpha, hatch=hatch, color=col, edgecolor=ec, fill=fill, label=fn)
        ax2.bar(x-baroffset, h, width=barwidth, alpha=alpha, hatch=hatch, color=col, edgecolor=ec, fill=fill, label=fn)

        if ylim:
            ylim = max(np.max(h[:i99]), ylim)
        else:
            ylim = np.max(h[:i99])

        xlim = max(xlim, i99)

        ax2.set_ylim([0, ylim])
        ax1.vlines(i05, 0, 1.0, color=col, edgecolor=ec, lw=3, linestyle=vls)
        ax2.vlines(i05, 0, 1.0, color=col, edgecolor=ec, lw=3, linestyle=vls)
        ax1.vlines(i95, 0, 1.0, color=col, edgecolor=ec, lw=3, linestyle=vls)
        ax2.vlines(i95, 0, 1.0, color=col, edgecolor=ec, lw=3, linestyle=vls)


    ax1.set_ylim([0, 1.0])
    ax1.set_xlim([0, xlim])
    ax1.set_ylabel("CDF[n_reads(DP)]")
    ax2.set_ylabel("H")
    ax2.set_xlabel("DP")
    ax1.legend(loc="lower right")

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
