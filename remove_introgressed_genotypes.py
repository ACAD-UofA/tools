#!/usr/bin/env python
# Copyright (c) 2015 Graham Gower <graham.gower@gmail.com>
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
"""
Modify introgressed genotypes to be unknown.
"""

from __future__ import print_function
import sys
import re
import collections

class ParseError(Exception):
    pass

def parse_hybrid_blocks(filename):
    """
    Accumulate hybrid intervals, grouped by sample.

    .xls has been converted to tab separated file.
    """

    hblocks = collections.defaultdict(list)
    lineno = 0

    with open(filename) as f:
        # skip header
        next(f)
        lineno += 1

        for line in f:
            lineno += 1
            line = line.rstrip("\r\n")
            fields = line.split("\t")

            #num_loci = fields[0]
            sample_id = fields[1]
            chrom = fields[2]
            #plot_id_new = fields[3]
            start = int(fields[4])
            end = int(fields[5])
            #copy = int(fields[6])

            # last six chars of the id match our ids
            sample = sample_id[-6:]

            try:
                chrom_int = int(chrom)
            except ValueError:
                raise RuntimeError("Unknown chromosome {}".format(chrom_str))

            if chrom_int > 22:
                # Eigensoft is set for humans, with chr 23/24 hard coded as X/Y,
                # chromosomes above 22 were translated in the snps file by adding 10.
                chrom_int += 10
            chrom = str(chrom_int)

            hblocks[sample].append((chrom, start, end))

    return hblocks

def parse_ind(filename):
    """
    Return list of samples.
    """
    re_s = re.compile(r"[ \t]+")
    lineno = 0
    samples = []

    with open(filename) as f:
        for line in f:
            lineno += 1
            fields = re_s.split(line.strip())

            if len(fields) != 3:
                raise ParseError("{}: line {}: expected 3 fields, got {}.".format(filename, lineno, len(fields)))

            sample_id = fields[0]
            #u = fields[1]
            group = fields[2]

            # take only the last six chars
            samples.append(sample_id[-6:])

    return samples

def yield_snp_loci(snp_filename):
    re_s = re.compile(r"[ \t]+")
    lineno = 0
    ref = alt = None
    with open(snp_filename) as f:
        for line in f:
            lineno += 1
            fields = re_s.split(line.strip())

            if len(fields) != 4 and len(fields) != 6:
                raise ParseError("{}: line {}: expected 4 or 6 fields, got {}.".format(snp_filename, lineno, len(fields)))

            snpid = fields[0]
            chrom = fields[1]
            gpos = fields[2]
            pos = int(fields[3])

            if len(fields) == 6:
                ref = fields[4]
                alt = fields[5]

            yield (snpid, chrom, gpos, pos, ref, alt)

def yield_geno_line(geno_filename):
    with open(geno_filename) as f:
        for line in f:
            line = line.rstrip("\r\n")
            yield line

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("usage: {} hybrid_blocks.tsv file.ind file.geno file.snp out.geno out.snp".format(sys.argv[0]))
        sys.exit(1)

    hb_filename = sys.argv[1]
    ind_filename = sys.argv[2]
    geno_filename = sys.argv[3]
    snp_filename = sys.argv[4]
    geno_out_filename = sys.argv[5]
    snp_out_filename = sys.argv[6]

    hblocks = parse_hybrid_blocks(hb_filename)
    samples = parse_ind(ind_filename)

    snps = yield_snp_loci(snp_filename)
    genotypes = yield_geno_line(geno_filename)

    introgressed_samples = set()

    with open(geno_out_filename, "w") as f_geno, \
         open(snp_out_filename, "w") as f_snp:
        for (snpid, chrom, gpos, pos, ref, alt), geno in zip(snps, genotypes):

            gts = []
            for c, sample in zip(geno, samples):
                if sample not in hblocks:
                    gts.append(c)
                    continue

                for block_chrom, start, end in hblocks[sample]:
                    if chrom == block_chrom and start <= pos <= end:
                        #print("hit: {}, {}, {}, {}, {}".format(sample, chrom, start, pos, end), file=sys.stderr)
                        # set genotype as unknown
                        gts.append("9")
                        introgressed_samples.add(sample)
                        break
                else:
                    gts.append(c)

            if sum(1 for gt in gts if gt!="9") == 0:
                # only missing genotypes remain
                continue

            if ref is not None:
                # 6 col
                print(snpid, chrom, gpos, pos, ref, alt, sep="\t", file=f_snp)
            else:
                # 4 col
                print(snpid, chrom, gpos, pos, sep="\t", file=f_snp)
            print("".join(gts), file=f_geno)


    #print("\n".join(introgressed_samples), file=sys.stderr)
