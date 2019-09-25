#!/usr/bin/env python
#
# Take a refseq *_genomic.fna.gz file and print the repeat mask as a bed file.
#
# From the refseq README.txt for any genome reference assembly:
# ``Repetitive sequences in eukaryotic genome assembly sequence files, as 
# identified by WindowMasker (Morgulis A, Gertz EM, Schaffer AA, Agarwala R. 
# 2006. Bioinformatics 22:134-41), have been masked to lower-case.''

import re

def parse_fa(filename):
    contig_header = None
    seq = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line[0] == ">":
                if contig_header:
                    yield contig_header, "".join(seq)
                contig_header = line
                seq = []
            else:
                seq.append(line)

    # last sequence
    yield contig_header, "".join(seq)


if __name__ == "__main__":
    re_space = re.compile(r"[ \t]+")
    for contig_header, seq in parse_fa("/dev/stdin"):
        contig = re_space.split(contig_header[1:])[0]
        start = end = None
        for i, s in enumerate(seq, 0):
            if s.islower():
                end = i+1
                if start is None:
                    start = i
            else:
                if start is not None:
                    print("{}\t{}\t{}".format(contig, start, end))
                    start = end = None
                    exit(1)
        if start is not None:
            print("{}\t{}\t{}".format(contig, start, end))
