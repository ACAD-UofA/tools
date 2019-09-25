#!/usr/bin/env python

from __future__ import print_function
import sys
from parsers import vcf_parser

def parse_geno(fn):
    with open(fn) as f:
        for line in f:
            yield line.strip()

def parse_ind(fn):
    with open(fn) as f:
        for line in f:
            fields = line.split()
            sample = fields[0]
            group = fields[2]
            yield sample, group

def parse_snp(fn):
    with open(fn) as f:
        for line in f:
            fields = line.split()
            #id = fields[0]
            chrom = fields[1]
            #gpos = fields[2]
            pos = int(fields[3])
            #ref = fields[4]
            #alt = fields[5]
            yield chrom, pos

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Convert VCF and EIGENSOFT files to SMC++ input format")
    parser.add_argument("--min-depth", type=int, default=5, help="sites with DP lower than this are marked unknown")
    parser.add_argument("--max-depth", type=int, default=1000, help="sites with DP greater than this are marked unknown")
    parser.add_argument("-g", "--group", type=str, help="restrict to GROUP in the EIGENSOFT files")
    parser.add_argument("vcf", help="input VCF for distinguished individual")
    parser.add_argument("eig_pfx", help="prefix of *.geno *.snp *.ind files")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    geno_fn = args.eig_pfx + ".geno"
    snp_fn = args.eig_pfx + ".snp"
    ind_fn = args.eig_pfx + ".ind"

    indlist = list(parse_ind(ind_fn))
    if args.group:
        ii = [i for i,(ind,gr) in enumerate(indlist) if gr==args.group]
        if len(ii) == 0:
            print("Error: no individuals in group ``{}''".format(args.group), file=sys.stderr)
            exit(1)
    else:
        ii = list(range(len(indlist)))

    alleles = {}
    for (chrom,pos),g in zip(parse_snp(snp_fn), parse_geno(geno_fn)):
        gts = [int(g[i]) for i in ii if g[i]!="9"]
        if len(gts) == 0:
            continue
        n_ref = sum(gts)
        n_alt = 2*len(gts) - n_ref
        alleles[(chrom,pos)] = (n_alt, n_ref+n_alt)

    vparser = vcf_parser(args.vcf)
    vcf_samples = next(vparser)

    assert len(vcf_samples)==1, "multi-sample vcf not supported"
    sample = vcf_samples[0]

    gt_map = {"./.":-1, ".|.":-1}
    for a0 in range(9):
        for a1 in range(9):
            n_alts = sum(1 if a else 0 for a in [a0, a1])
            gt_map["{}/{}".format(a0,a1)] = n_alts
            gt_map["{}|{}".format(a0,a1)] = n_alts

    prev_pos = 0
    prev_chrom = None
    last = (-1, -1, -1)
    span = 0
    for row in vparser:
        lineno, _, chrom, pos, _, ref, alt, _, _, _, _, genotypes = row
        if len(ref) != len(alt):
            # skip indels, etc.
            continue

        gt_str = genotypes[sample]["GT"]
        try:
            d = gt_map[gt_str]
        except IndexError:
            print("{}: {}: unexpected genotype ``{}''".format(args.vcf, lineno, gt_str), file=sys.stderr)
            d = -1

        try:
            dp = int(genotypes[sample]["DP"])
        except (IndexError, ValueError):
            print("{}: {}: missing FORMAT/DP field".format(args.vcf, lineno), file=sys.stderr)
            exit(1)

        if dp < args.min_depth or dp > args.max_depth:
            d = -1

        if chrom.startswith("chr"):
            chrom = chrom[3:]

        u1, n1 = alleles.get((chrom,pos), (0, 0))

        #print(span, d, u1, n1, sep="\t", file=sys.stderr)

        if prev_chrom != chrom:
            if prev_chrom is None:
                prev_crhom = chrom
            else:
                print("{}: {}: vcf contains more than one chromosome, bailing".format(args.vcf, lineno), file=sys.stderr)
                exit(1)

        if prev_pos+1 != pos:
            # fill gap to prev_pos
            if span > 0 and last != (-1, 0, 0):
                print(span, *last, sep="\t")
                span = 0
            span += pos - prev_pos - 1
            last = (-1, 0, 0)

        if (d, u1, n1) == last:
            span += 1
        else:
            if span > 0:
                print(span, *last, sep="\t")
            last = (d, u1, n1)
            span = 1

        prev_pos = pos

    # mop up
    if span > 0:
        print(span, *last, sep="\t")
