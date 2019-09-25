
from __future__ import print_function
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import numpy as np
from numpy import random
import sys
import collections

from seqcontent import entropy, selfanneal, seq2probes

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


def parse_seqcontent_bed(filename, sc=None):
    if sc is None:
        sc = collections.defaultdict(list)
    with open(filename) as f:
        header = next(f)
        hfields = header.strip().split("\t")
        for line in f:
            fields = line.rstrip().split("\t")
            if len(hfields) != len(fields):
                continue
            for h, f in zip(hfields, fields):
                sc[h].append(f)
    return hfields, sc

def randseq(size):
    seq = random.choice(("A","C","G","T"), size=size)
    return "".join(seq)

def randalt(x):
    if x is "A":
        return random.choice(("C","G","T"))
    if x is "C":
        return random.choice(("A","G","T"))
    if x is "G":
        return random.choice(("A","C","T"))
    if x is "T":
        return random.choice(("A","C","G"))

def snp_stats(sc):
    snp_counts = collections.Counter(zip(sc["ref"], sc["alt"]))
    tiset = set([("C","T"), ("T","C"), ("G","A"), ("A","G")])
    ti = 0
    tv = 0
    print("ref\talt\tcount")
    for ref in "ACGT":
        for alt in "ACGT":
            if ref==alt:
                continue

            n = snp_counts[(ref,alt)]
            if (ref,alt) in tiset:
                ti += n
            else:
                tv += n

            print(ref, alt, n, sep="\t")

    print("\n\n")
    print("Ti", ti, sep="\t")
    print("Tv", tv, sep="\t")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("usage: {} file1.bed [... fileN.bed]".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    sc = None
    for fn in sys.argv[1:]:
        hfields, sc = parse_seqcontent_bed(fn, sc)

    snp_stats(sc)

    plen = 60

    # simulate nreps sequences
    nreps = 1000
    e = np.empty(nreps)
    sa = np.empty(nreps)
    for i in range(nreps):
        s = randseq(121)
        ref = s[plen+1]
        alt = randalt(ref)
        e_min = float('inf')
        sa_max = 0
        for p in seq2probes(s, ref, alt, plen):
            e_min = min(entropy(p), e_min)
            sa_max = max(selfanneal(p), sa_max)
        e[i] = e_min
        sa[i] = sa_max


    setCustomHatchWidth(1)

    e_probes = map(float, sc["entropy"])
    sa_probes = map(float, sc["selfanneal"])

    pdf = PdfPages("seqcontent.pdf")
    #fig_w, fig_h = plt.figaspect(3.0/4.0)
    #fig1 = plt.figure(figsize=(fig_w, fig_h))
    #fig1 = plt.figure(figsize=(8.27,11.69)) # a4 portrait
    fig1 = plt.figure(figsize=(11.69,8.27)) # a4 landscape
    gs1 = gridspec.GridSpec(2, 1)
    ax1, ax2 = [fig1.add_subplot(gs) for gs in gs1]

    e_bound1 = np.percentile(e_probes, 5)
    e_bound2 = np.percentile(e, 5)
    sa_bound1 = np.percentile(sa_probes, 95)
    sa_bound2 = np.percentile(sa, 95)

    bins1 = np.arange(2.5,4.0,(4.0-2.5)/50)
    ax1.hist(e_probes, bins=bins1, normed=True, label='probes', color="lightblue", edgecolor='darkblue', hatch='', fill=True)
    ax1.hist(e, bins=bins1, normed=True, label='random', color="#660000", edgecolor='#660000', hatch='xx', fill=False)
    ylim = ax1.get_ylim()
    ax1.vlines(e_bound1, 0, ylim[1], color="lightblue", edgecolor='darkblue', lw=3, linestyle='--')
    ax1.vlines(e_bound2, 0, ylim[1], color='#660000', edgecolor='#660000', lw=3, linestyle='--')
    ax1.set_title('Shannon\'s Entropy (3-mer symbols)')
    ax1.set_ylabel('density')
    ax1.set_xlabel('$H_3$')
    ax1.legend(loc='upper left')
    
    bins2 = list(range(35))
    ax2.hist(sa_probes, bins=bins2, normed=True, label='probes', color="lightblue", edgecolor='darkblue', hatch='', fill=True)
    ax2.hist(sa, bins=bins2, normed=True, label='random', color="#660000", edgecolor='#660000', hatch='xx', fill=False)
    ylim = ax2.get_ylim()
    ax2.vlines(sa_bound1, 0, ylim[1], color="lightblue", edgecolor='darkblue', lw=3, linestyle='--')
    ax2.vlines(sa_bound2, 0, ylim[1], color='#660000', edgecolor='#660000', lw=3, linestyle='--')
    ax2.set_title('Self annealing property')
    ax2.set_ylabel('density')
    ax2.set_xlabel('Smith-Watterman alignment score')
    ax2.legend()

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
