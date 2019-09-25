
from __future__ import print_function
import sys
import collections
import math
import swalign
import numpy as np

from parsers import parse_fa, vcf_parser

def kmers(seq, k):
    """
    yield all kmers in seq
    """
    for i in range(len(seq)-k+1):
        kmer = seq[i:i+k]
        yield kmer

def entropy(seq, k=3):
    """
    Shannon's entropy, where symbols are kmers in the sequence.
    """
    counts = collections.Counter(kmers(seq, k))
    l = np.sum(counts.values())
    h = 0
    for x, n in counts.iteritems():
        px = float(n)/l
        h += px * math.log(px)
    return -h

scoring = swalign.NucleotideScoringMatrix()
sw = swalign.LocalAlignment(scoring)

def selfanneal(seq):
    """
    Smith-Waterman local alignment score of seq to its reverse complement
    """
    aln = sw.align(seq, swalign.revcomp(seq))
    return aln.score

def gc(seq):
    """
    GC fraction
    """
    counts = collections.Counter(seq)
    return float(counts["C"]+counts["G"])/np.sum(counts.values())

def cpg(seq):
    """
    Number of CpGs
    """
    l = len(seq)
    if l < 2:
        return 0
    n = 0
    for i in range(l-1):
        if seq[i:i+2] == "CG":
            n += 1
    return n

def windowmasker(seq):
    """
    Number of sites masked (lowercase, or non ACGT) in the seq
    """
    counts = collections.Counter(seq)
    unmasked = np.sum(counts[x] for x in "ACGT")
    return len(seq) - unmasked

def seqcontent(seq, ref, alt, plen):
    """
    Obtain some metrics about the seq
    """
    probes = seq2probes(seq, ref, alt, plen)
    h_a = []
    gc_a = []
    sa_a = []
    cpg_a = []
    for probe in probes:
        h_a.append(entropy(probe))
        gc_a.append(gc(probe))
        sa_a.append(selfanneal(probe))
        cpg_a.append(cpg(probe))

    return np.min(h_a), np.mean(gc_a), np.max(sa_a), np.mean(cpg_a)

def seq2probes(seq, ref, alt, plen):
    """
    Return a list of tiled probes for the seq
    """
    assert(len(seq) == plen*2 +1)
    assert(ref != alt)
    seq = seq.upper()
    plen2 = plen/2
    probes = [seq[:plen],
                seq[plen2:plen+plen2+1],
                seq[plen2:plen] + alt + seq[plen+1:plen+plen2+1],
                seq[plen+1:]]
    return probes

def parse_mappability(filename):
    with open(filename) as f:
        length = int(next(f))
        mapmap = np.array(length)
        i = 0
        for line in f:
            line = line.rstrip("\n")
            n, x = line.split(" ")
            mapmap[i:i+n] = x
            i += n
    return mapmap

def parse_faidx(filename):
    offset = 0
    faidict = {}
    with open(filename) as f:
        for line in f:
            line = line.rstrip("\n")
            fields = line.split("\t")
            chrom = fields[0]
            length = int(fields[1])
            faidict[chrom] = offset
            offset += length 

_mapmap = None
_faidict = None
def mappability(chrom, start, end, mapfile, faidx, k=23):
    """
    Return mean mappability over the region described by chrom:start-end.
    """
    global _mapmap
    global _faidict
    if _mapmap is None:
        _mapmap = parse_mappability(mapfile)
    if _faidict is None:
        _faidict = parse_faidx(faidx)

    ch_offset = _faidict[chrom]
    m = [_mapmap[x] for x in range(ch_offset+start, ch_offset+end+1-k)]
    return np.mean(m)

def filter_snps_by_dist(vcf_filename, plen):
    """
    Yield SNPs which are plen positions away from other SNPs.
    """

    vp = vcf_parser(vcf_filename, yield_samples=True, parse_genotypes=True)
    samples = next(vp)
    if len(samples) != 1:
        raise Exception("{} samples. Single sample vcf expected.".format(len(samples)))

    sample = samples[0]

    last  = None
    exclude = False

    for vline in vp:
        chrom = vline[2]
        pos = vline[3]
        ref = vline[5].upper()
        alt = vline[6].upper()
        genotype_field = vline[11][sample]
        ad = genotype_field["AD"]
        dp = np.sum(map(int, ad.split(",")))

        if last is not None:
            if last[0] == chrom and pos-last[1] < plen:
                exclude = True
            elif not exclude:
                yield last
            else:
                exclude = False

        last = [chrom, pos, ref, alt, dp]

    if not exclude:
        yield last
            

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: {} in.fa in.vcf".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    fp = parse_fa(sys.argv[1])

    snpdict = {}

    for label, seq in fp:
        chrom, extent = label[1:].split(":")
        start, end = map(int, extent.split("-"))
        pos = start + (end-start) / 2
        if windowmasker(seq) > 0:
            continue
        snpdict[(chrom, pos)] = seq

    plen = 60

    print("chrom", "pos", "ref", "alt",
            "downstream60bp", "upstream60bp", "DP",
            "entropy", "GC", "selfanneal", "CpGs",
            sep="\t")

    for chrom, pos, ref, alt, dp in filter_snps_by_dist(sys.argv[2], plen):
        if len(alt) != 1:
            continue
        if (chrom, pos) in snpdict:
            seq = snpdict[(chrom, pos)]
            sc = seqcontent(seq, ref, alt, plen)
            print(chrom, pos, ref, alt, seq[:plen], seq[-plen:], dp, *sc, sep="\t")
