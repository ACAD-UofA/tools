#!/usr/bin/env python
"""
Append gene annotation columns to a tsv file.
"""

from __future__ import print_function
import sys
import collections
from intervaltree import IntervalTree

def parse_annot(filename):
    lineno = 0

    ensembl_id_annot = {}
    intervals = collections.defaultdict(IntervalTree)

    with open(filename) as f:

        # get headers
        line = next(f).rstrip("\r\n")
        headers = line.split("\t")

        for line in f:
            line = line.rstrip("\r\n")
            fields = line.split("\t")

            if len(headers) != len(fields):
                raise ParseError("{}: line {}: expected {} columns, but got {}.".format(filename, lineno, len(headers), len(fields)))

            gene = {}
            for k, v in zip(headers, fields):
                gene[k] = v

            chrom = gene["chrom"]
            start = int(gene["start"])
            end = int(gene["end"])
            ensembl_id = gene["ensembl_id"]

            if chrom.startswith("chr"):
                chrom = chrom[3:]

            ensembl_id_annot[ensembl_id] = gene
            intervals[chrom][start:end] = ensembl_id

    return headers, intervals, ensembl_id_annot

def parse_tsv(filename):
    with open(filename) as f:

        # header
        line = next(f)
        yield line.rstrip("\r\n")

        data = []

        for line in f:
            line = line.rstrip("\r\n")
            fields = line.split("\t")
            chrom = fields[0]
            pos = int(fields[1])

            if chrom.startswith("chr"):
                chrom = chrom[3:]

            yield (chrom, pos, line)


def parse_args():
    import argparse

    parser = argparse.ArgumentParser(description="Append gene annotation columns to a tsv file.")
    parser.add_argument("tsv_file", help="input tsv with CHROM, POS as first two columns")

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    headers, annot_intervals, ensembl_id_annot = parse_annot("/localscratch/grg/H_sapiens/data/annot.tsv")

    headers = headers[1:] # dont add the chromosome, its already there

    tparser = parse_tsv(args.tsv_file)

    tsv_header = next(tparser)

    sys.stdout.write(tsv_header)
    for header in headers:
        sys.stdout.write("\t{}".format(header))
    sys.stdout.write("\n")

    for chrom, pos, line in tparser:
        sys.stdout.write(line)

        intervals = annot_intervals[chrom][pos]

        if len(intervals) >= 1:
            #if len(intervals) > 1:
            #    print("Warning: multiple genes for {}:{}".format(chrom, pos), file=sys.stderr)
            ensembl_id = list(intervals)[0].data
            annot = ensembl_id_annot[ensembl_id]

            for header in headers:
                sys.stdout.write("\t{}".format(annot[header]))
        else:
            sys.stdout.write("\t"*len(headers))

        sys.stdout.write("\n")
