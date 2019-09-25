#!/usr/bin/env python

# https://www.cog-genomics.org/plink2/formats

from __future__ import print_function
import sys
import itertools
import collections
from scipy.stats import beta

import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

from colours import pal16

class InvalidBed(Exception):
    pass
class InvalidLocus(Exception):
    pass

def parse_fam(fn):
    with open(fn) as f:
        for line in f:
            fields = line.split()
            fid = fields[0] # family id
            iid = fields[1] # individual id
            #pid = fields[2] # paternal id (0 otherwise)
            #mid = fields[3] # maternal id (0 otherwise)
            #sex = fields[4] # 1=male, 2=female, 0=unknown
            #pheno = fields[5] # e.g. 1=control, 2=case, 0=missing
            yield fid, iid

def parse_bim(fn):
    with open(fn) as f:
        for line in f:
            fields = line.split()
            chrom = fields[0]
            #rs = fields[1]
            #gpos = fields[2] # genetic position in centimorgans (0 otherwise)
            pos = int(fields[3])
            #a1 = fields[4]
            #a2 = fields[5]
            yield chrom, pos

def parse_bed(fn, bufsize=16*1024):
    with open(fn, "rb") as f:
        magic = f.read(3)
        if magic != b"\x6c\x1b\x01":
            x = [hex(ord(m)) for m in magic]
            raise InvalidBed("{}: invalid magic ({},{},{})".format(fn,*x))

        while True:
            buf = f.read(bufsize)
            if not buf:
                break
            for c in buf:
                b = ord(c)
                for bit in range(4):
                    yield (b>>bit)&0x3

class PlinkBfile():
    def __init__(self, bed_fn, bim_fn=None, fam_fn=None, magic=b"\x6c\x1b\x01"):

        if not bed_fn.endswith(".bed"):
            self.bed_fn = bed_fn + ".bed"
            if bim_fn is None:
                self.bim_fn = bed_fn + ".bim"
            else:
                self.bim_fn = bim_fn
            if fam_fn is None:
                self.fam_fn = bed_fn + ".fam"
            else:
                self.fam_fn = fam_fn
        else:
            self.bed_fn = bed_fn

        with open(self.bed_fn, "rb") as f:
            self.msize = len(magic)
            m = f.read(self.msize)
            if m != magic:
                x = [hex(ord(mi)) for mi in m]
                raise InvalidBed("{}: invalid magic ({},{},{})".format(self.bed_fn,*x))

        self.indlist = list(parse_fam(self.fam_fn))
        self.loci = None

    def __enter__(self):
        return self._gen_dataframe()

    def _gen_dataframe(self):
        locusparser = parse_bim(self.bim_fn)
        gtparser = parse_bed(self.bed_fn)
        n = len(self.indlist)
        for chrom, pos in locusparser:
            gts = itertools.islice(gtparser, n)
            for (fid, iid), gt in zip(self.indlist, gts):
                yield chrom, pos, fid, iid, gt

    def _get_loci(self):
        if self.loci is not None:
            return
        self.loci = {(ch,pos):i for i,(ch,pos) in enumerate(parse_bim(self.bim_fn))}

    def locus_gts(self, chrom, pos):
        self._get_loci()
        i = self.loci.get((chrom,pos), None)
        if i is None:
            raise InvalidLocus("couldn't find {}:{}".format(chrom,pos))

        n = len(self.indlist)
        blksize = n/4
        if n%4 != 0:
            blksize += 1

        with open(self.bed_fn) as f:
            f.seek(self.msize + i*blksize)
            buf = f.read(blksize)
            if len(buf) != blksize:
                raise InvalidLocus("expected {} bytes, but only got {}".format(blksize, len(buf)))
            j = 0
            for c in buf:
                b = ord(c)
                for bit in range(4):
                    j += 1
                    if (j > n):
                        return
                    yield (b>>bit)&0x3


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: {} bfile chrI:posI [...chrN:posN]".format(sys.argv[0]),
                file=sys.stderr)
        exit(1)

    bfile = sys.argv[1]

    pb = PlinkBfile(bfile)

    families = collections.defaultdict(list)
    for i, (fid, _) in enumerate(pb.indlist):
        families[fid].append(i)

    gtmap = {0:0, 1:-1, 2:1, 3:2}
    #print("FID", "AC", "N", "AF", "95%CI", sep="\t")

    pdf = PdfPages("out.pdf")
    fig_w, fig_h = plt.figaspect(9.0/16.0)
    fig1 = plt.figure(figsize=(fig_w,fig_h))
    gs1 = gridspec.GridSpec(1, 1)
    ax1 = fig1.add_subplot(gs1[0])

    x = list(range(1,len(families.keys())+1))

    for arg in sys.argv[2:]:
        chrom, pos = arg.split(":")
        gts = list(pb.locus_gts(chrom, int(pos)))
        aflist = []
        ci1list = []
        ci2list = []
        for fid, ilist in families.iteritems():
            gt = filter(lambda g: g!=-1, (gtmap[gts[i]] for i in ilist))
            ac = sum(gt)
            n = 2*len(gt) # some may be haploidised, so this is skewed
            if n == 0:
                af = float('nan')
            else:
                af = float(ac)/n

            # Beta 95%CI with non-informative prior, aka Jefferey's interval.
            # See Brown, Cai, and DasGupta (2001). doi:10.1214/ss/1009213286
            ci = beta.interval(0.95, ac+0.5, (n-ac)+0.5)

            #print(fid, ac, n, af, ci, sep="\t")

            aflist.append(af)
            # this is dumb
            ci1list.append(af-ci[0])
            ci2list.append(ci[1]-af)

        ax1.errorbar(x, aflist, yerr=[ci1list,ci2list], color="black", lw=1.5, capsize=0, fmt="o")

        ax1.set_xticks(x)
        ax1.set_xticklabels(families.keys(), rotation="45", ha="right", fontsize=10)
        ax1.set_xlim([min(x)-0.5, max(x)+0.5])
        ax1.set_ylabel("Allele Frequency")

        title = arg
        if not title.startswith("chr"):
            title = "chr" + title
        ax1.set_title(title)

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
