#!/usr/bin/env python

from __future__ import print_function
import sys
import os.path
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import numpy as np


def setCustomHatchWidth(customWidth):
    """
    Monkeypatch the pdf writing for hatches, to change the hatch line width
    """

    # make sure you have the correct imports,
    # they may differ depending on the matplotlib version
    import matplotlib.backends.backend_pdf
    #from matplotlib.externals import six
    from matplotlib.backends.backend_pdf import Name, Op
    from matplotlib.transforms import Affine2D
    from matplotlib.path import Path

    def _writeHatches(self):
        hatchDict = dict()
        sidelen = 72.0
        for hatch_style, name in self.hatchPatterns.iteritems():
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


def parse_input(fn, args):
    data = []
    skip = args.skip
    with open(fn) as f:
        while skip > 0:
            next(f)
            skip -= 1
        for line in f:
            fields = line.split()
            if len(fields) < args.column:
                continue
            val = float(fields[args.column-1])
            if np.isinf(val) or np.isnan(val):
                continue
            data.append(float(fields[args.column-1]))

    if args.integer:
        return np.array(data, dtype=int)
    elif args.float:
        return np.array(data, dtype=float)

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="plot histogram for 1 column of input file")
    parser.add_argument("-c", "--column", type=int, required=True, help="1-based column to plot")
    parser.add_argument("-s", "--skip", type=int, default=0, help="number of input lines to skip [%(default)s]")
    parser.add_argument("-b", "--bins", type=int, help="number of histogram bins [%(default)s]")
    parser.add_argument("-i", "--integer", action="store_true", default=False, help="column is an integer [%(default)s]")
    parser.add_argument("-f", "--float", action="store_true", default=True, help="column is a float [%(default)s]")
    parser.add_argument("-n", "--normalise", action="store_true", default=False, help="normalise histogram so area=1 [%(default)s]")
    parser.add_argument("-w", "--wide", action="store_true", default=False, help="plot widescreen ratio (16x9) [%(default)s]")
    parser.add_argument("-l", "--logy", action="store_true", default=False, help="plot y-axis on log scale [%(default)s]")
    parser.add_argument("-o", "--opdf", type=str, default="out.pdf", help="output filename [%(default)s]")
    parser.add_argument("-t", "--title", type=str, help="plot title")
    parser.add_argument("-x", "--xlabel", type=str, help="x axis label")
    parser.add_argument("-y", "--ylabel", type=str, help="y axis label")
    parser.add_argument("infiles", nargs="+", help="input file(s)")
    args = parser.parse_args()

    if args.integer:
        args.float = False

    return args


if __name__ == "__main__":

    args = parse_args()

    setCustomHatchWidth(1)

    data = {}
    for fn in args.infiles:
        data[fn] = parse_input(fn, args)

    pdf = PdfPages(args.opdf)
    if args.wide:
        fig_w, fig_h = plt.figaspect(9.0/16.0)
    else:
        fig_w, fig_h = plt.figaspect(3.0/4.0)
    fig1 = plt.figure(figsize=(fig_w, fig_h))
    gs1 = gridspec.GridSpec(1, 1)
    ax1 = fig1.add_subplot(gs1[0])

    if not args.title:
        if len(args.infiles) == 1:
            args.title = os.path.basename(args.infiles[0])
        else:
            args.title = ""

    if not args.ylabel:
        if args.normalise:
            if args.logy:
                args.ylabel = "log(Density)"
            else:
                args.ylabel = "Density"
        else:
            if args.logy:
                args.ylabel = "log(Counts)"
            else:
                args.ylabel = "Counts"
    if not args.xlabel:
        args.xlabel = "column {}".format(args.column)

    if not args.bins:
        if args.integer:
            args.bins = np.arange(np.min(d[fn] for fn in data), np.max(data)+1)
        else:
            args.bins = 30


    pal = ["#006050", "#806070", "#959025"]
    facecols = ['lightblue', '#660000', 'yellow']
    edgecols = ['darkblue', '#660000', 'darkyellow']
    hatches = ['', 'xx', '--']
    fills = [True, False, True]

    for fn, fcol, ecol, hat, fill in zip(args.infiles, facecols, edgecols, hatches, fills):
        if args.bins:
            bins = args.bins
        else:
            if args.integer:
                bins = np.arange(np.min(data[fn]), np.max(data[fn])+1)
            else:
                bins = None

        #bins = np.arange(-0.5, 0.5, 1.0/100)
        ax1.hist(data[fn], bins=bins, normed=args.normalise, color=fcol, edgecolor=ecol, hatch=hat, fill=fill, label=os.path.basename(fn))

    ax1.set_title(args.title)
    ax1.set_xlabel(args.xlabel)
    ax1.set_ylabel(args.ylabel)

    if len(args.infiles) > 1:
        ax1.legend(loc='upper left')

    if args.logy:
        plt.yscale('log', nonposy='clip')



    plt.tight_layout()
    pdf.savefig()
    pdf.close()
