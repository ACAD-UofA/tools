#!/usr/bin/env python

from __future__ import print_function
import sys
from itertools import ifilter

# https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
def parse_gff(filename):
    if filename.endswith(".gz"):
        import gzip
        xopen = gzip.open
    else:
        xopen = open

    with xopen(filename) as f:
        for line in f:
            if line[0] == "#":
                continue
            line = line.rstrip()
            fields = line.split("\t")
            seqid = fields[0]
            #source = fields[1]
            type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            #score = fields[5]
            #strand = fields[6]
            #phase = fields[7]
            attr = fields[8]

            yield seqid, type, start, end, attr

def parse_assembly_report(filename):
    with open(filename) as f:
        for line in f:
            if line[0] == "#":
                continue
            line = line.rstrip()
            fields = line.split("\t")
            #seqname = fields[0]
            #seqrole = fields[1]
            molecule = fields[2]
            #molecule_type = fields[3]
            genbank = fields[4]
            #relationship = fields[5]
            refseq = fields[6]
            #assembly = fields[7]
            #seqlen = fields[8]
            #ucsc = fields[9]
            yield molecule, genbank, refseq

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="convert GFF3 to BED")
    parser.add_argument("-t", "--type", action="append", help="extract this 'type' of entry from the GFF")
    parser.add_argument("gff", help="genomic.gff.gz input file")
    parser.add_argument("report", help="assembly_report.txt")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()

    refmap = {}
    for molecule, genbank, refseq in parse_assembly_report(args.report):
        if molecule != "na":
            refmap[refseq] = molecule
        elif genbank != "na":
            refmap[refseq] = genbank

    if args.type is None:
        gffparser = parse_gff(args.gff)
    else:
        types = set(args.type)
        predicate = lambda x: x[1] in types
        gffparser = ifilter(predicate, parse_gff(args.gff))

    for seqid, _, start, end, _ in gffparser:
        chrom = refmap[seqid]
        print(chrom, start-1, end, sep="\t")
