
from __future__ import print_function
import sys
import collections
import math
import numpy as np

def kmers(seq, k):
    for i in range(len(seq)-k+1):
        kmer = seq[i:i+k]
        if "N" not in kmer:
            yield kmer

def centropy(seq, counts=None):
    """
    Conditional entropy H(X|Y) for nucleotide X given previous nucleotide Y.
    """
    if counts is None:
        counts = collections.Counter(seq)
    dicounts = collections.Counter(kmers(seq, 2))
    l = len(seq)
    dil = np.sum(dicounts.values())
    h = 0
    for y in "ACGT":
        ny = counts[y]
        if ny == 0:
            continue
        py = float(ny)/l        # p(y)
        for x in "ACGT":
            nxy = dicounts[y+x]
            if nxy == 0:
                continue
            pxy = float(nxy)/dil    # p(x|y)
            h += py * pxy * math.log(pxy)
    return -h

def entropy(seq, k=1):
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

from numpy import random
def randseq(size):
    seq = random.choice(("A","C","G","T"), size=size)
    return "".join(seq)

if __name__ == "__main__":
    for s in sys.argv[1:]:
        s = s.upper()
        print(s, centropy(s),entropy(s,3),entropy(s,4), sep="\t")


    c = []
    e = []
    z = []

    for x in range(1000):
        s = randseq(60)
        c.append(centropy(s))
        e.append(entropy(s,3))
        z.append(entropy(s,4))

    print("mean", np.mean(c), np.mean(e), np.mean(z), sep="\t")
    print("var", np.var(c), np.var(e), np.var(z), sep="\t")
    print("mean-2SD", np.mean(c)-2*np.std(c), np.mean(e)-2*np.std(e), np.mean(z)-2*np.std(z), sep="\t")
