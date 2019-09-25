#!/usr/bin/env python

from __future__ import print_function
import sys
import numpy as np

def parse_dss(fn):
    if fn.endswith(".gz"):
        import gzip
        xopen = gzip.open
    else:
        xopen = open

    with xopen(fn) as f:
        for line in f:
            fields = line.split("\t")

            chrom = fields[0]
            if not chrom.startswith("chr"):
                # ignore scaffolds
                continue
            chrom = chrom[len("chr"):]
            try:
                chrom = int(chrom)
            except ValueError:
                # ignore sex chromosomes
                continue

            pos = int(fields[1])
            n = int(fields[2])
            x = int(fields[3])

            yield chrom, pos, float(x)/n

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: {} dss.txt".format(sys.argv[0]), file=sys.stderr)
        sys.exit(1)

    dss_fn = sys.argv[1]
    dss_parser = parse_dss(dss_fn)
    maxdist = 10000
    binsize = 10

    mean = 0.0
    n = 0
    var_ss = 0.0
    cov_ss = np.zeros(maxdist/binsize, dtype=float)
    cov_ss_n = np.zeros(maxdist/binsize, dtype=float)

    chrom, pos, m = next(dss_parser)
    chrom_v = [chrom,]
    pos_v = [pos,]
    m_v = [m,]

    while len(pos_v) > 0:
        chrom0 = chrom_v.pop(0)
        pos0 = pos_v.pop(0)
        m0 = m_v.pop(0)

        n += 1
        last_mean = mean
        mean += (m0-mean)/n
        var_ss += (m0-mean)*(m0-last_mean)

        while True:
            try:
                chrom, pos, m = next(dss_parser)
            except StopIteration:
                break
            chrom_v.append(chrom)
            pos_v.append(pos)
            m_v.append(m)
            if chrom != chrom0 or pos-pos0 >= maxdist:
                break

        for chrom, pos, m in zip(chrom_v, pos_v, m_v):
            if chrom != chrom0 or pos-pos0 >= maxdist:
                break
            d = pos-pos0
            d /= binsize
            cov_ss_n[d] += 1
            # Bennett et al. 2009,  DOI: 10.1109/CLUSTR.2009.5289161
            # Eq III.9
            cov_ss[d] += (cov_ss_n[d]-1) * (m-mean)*(m0-mean) / cov_ss_n[d]

    var = var_ss / (n-1)
    #print(n, mean, var, sep="\t")

    print("delta", "n", "cov", "cor", sep="\t")
    for i in range(1,maxdist/binsize):
        if cov_ss_n[i] > 100:
            cov = cov_ss[i] / (cov_ss_n[i]-1)
            p = cov / var
        else:
            break
        print(i*binsize, int(cov_ss_n[i]), cov, p, sep="\t")
