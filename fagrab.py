#!/usr/bin/env python

from __future__ import print_function
import sys

def print_seq(contig_header, seq, linelen=70):
    print(contig_header)
    i = 0
    while i<len(seq):
        print(seq[i:i+linelen])
        i += linelen

class ParseError(Exception):
    pass

def parse_region(region):
    rfields = region.split(":")
    if len(rfields) > 2:
        raise ParseError("too many colons `{}'".format(region))

    chrom = rfields[0]
    start = 0
    end = -1
    if len(rfields) == 2:
        pfields = rfields[1].split("-")
        if len(pfields) != 2:
            raise ParseError("too many hyphens `{}'".format(region))
        try:
            start, end = int(pfields[0]), int(pfields[1])
        except ValueError:
            raise ParseError("invalid region, expected integers `{}'".format(rfields[1]))

    return chrom, start, end

def fasta_extract(fn, chrom, start, end):

    in_chr = True
    slist = []

    with open(fn) as f:
        for line in f:
            if line[0] == ">":
                if line[1:].startswith(chrom):
                    in_chr = True
                    continue
                elif in_chr:
                    # we moved out of the requested chromosome
                    break
            if in_chr:
                line = line.rstrip()
                slist.append(line.rstrip())

    
    seq = "".join(slist)
    if end == -1:
        return seq
    else:
        return seq[start:end]

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: {} in.fasta region".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    fasta_fn = sys.argv[1]
    region = sys.argv[2]

    chrom, start, end = parse_region(region)

    seq = fasta_extract(fasta_fn, chrom, start, end)
    print_seq(">{}".format(region), seq)
