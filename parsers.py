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


from __future__ import print_function
import sys
from operator import itemgetter
import collections
from gzopen import gzopen

class ParseError(Exception):
    pass

def parse_fa(filename):
    """
    fasta or fastq
    """

    fasta = None
    lineno = 0
    state = 0
    label = None
    qual = None

    with gzopen(filename) as f:
        for line in f:
            lineno += 1
            line = line.rstrip("\r\n")
            if not line:
                continue

            if fasta is None:
                # determine file type
                if line[0] == ">":
                    fasta = True
                elif line[0] == "@":
                    fasta = False
                else:
                    raise ParseError("{}: line {}: did not recognise fasta nor fastq".format(filename, lineno))

            if fasta:
                if line[0] == ">":
                    if label:
                        yield label, "".join(seq)
                    label = line
                    seq = []
                else:
                    seq.append(line)
            else:
                # fastq
                if line[0] == "@":
                    if label:
                        yield label, "".join(seq), "".join(qual)
                    state = 1
                    label = line
                    seq = []
                    qual = []
                    continue
                elif line[0] == "+":
                    state = 2
                    continue

                if state == 1:
                    seq.append(line)
                elif state == 2:
                    qual.append(line)

    # last one
    if fasta:
        yield label, "".join(seq)
    else:
        yield label, "".join(seq), "".join(qual)

def vcf_parser(filename, yield_headers=False, yield_samples=True, yield_genotypes=True, parse_genotypes=True):
    """
    Parse vcf. Pyvcf does not meet my needs.  Yields specified vcf lines.

    Assumes the format is the same for every line.
    """

    samples = None

    with gzopen(filename) as f:

        # parse header
        for lineno, line in enumerate(f, 1):
            if yield_headers:
                yield line.rstrip("\r\n")
            if not line.startswith("#"):
                break
            if line.startswith("#CHROM"):
                # column header line
                # #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample0 [sample1 ...]
                line = line.rstrip("\r\n")
                fields = line.split("\t")
                if len(fields) < 9:
                    raise ParseError("{}: line {}: expected at least 9 columns for vcf column header.".format(filename, lineno))
                samples = fields[9:]
                break

        if yield_samples:
            if samples is None:
                raise ParseError("{}: no vcf column header found.".format(filename))
            yield samples

        if not yield_genotypes:
            return

        genotypes = collections.defaultdict(dict)
        fmt_fields = None

        # parse variants
        for line in f:
            lineno += 1
            line = line.rstrip("\r\n")
            fields = line.split("\t")

            chrom = fields[0]
            pos = int(fields[1])
            id = fields[2]
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            filter = fields[6]
            info_str = fields[7]
            format_str = fields[8]
            genotype_str_list = fields[9:]

            if chrom.startswith("chr"):
                chrom = chrom[3:]

            if not parse_genotypes:
                yield (lineno, chrom, pos, id, ref, alt, qual, filter, info_str, format_str, genotype_str_list)
                continue

            if fmt_fields is None:
                fmt_fields = format_str.split(":")

            for sample, genotype_str in zip(samples, genotype_str_list):
                gen_fields = genotype_str.split(":")
                for fmt, gen in zip(fmt_fields, gen_fields):
                    genotypes[sample][fmt] = gen

            yield (lineno, line, chrom, pos, id, ref, alt, qual, filter, info_str, fmt_fields, genotypes)


def parse_bed(filename, sort_by_id=False, sort_by_coord=False):

    loci_set = set()
    loci_list = []
    lineno = 0

    with open(filename) as f:
        for line in f:
            lineno += 1
            line = line.rstrip("\r\n")
            fields = line.split("\t")

            chrom = fields[0]
            #start = int(fields[1])
            end = int(fields[2])
            snpid = fields[3]

            if chrom.startswith("chr"):
                chrom = chrom[3:]

            if (chrom, end) in loci_set:
                print("Warning: {}: line {}: duplicate SNP id {} for locus {}:{}, ignoring.".format(filename, lineno, snpid, chrom, end), file=sys.stderr)
                continue

            loci_set.add((chrom, end))
            loci_list.append((chrom, end, snpid))

    if sort_by_id:
        # Sort by SNP id.
        loci_list.sort(key=itemgetter(2))
    elif sort_by_coord:
        # sort by genomic coordinate
        def chr_cmp(a, b):
            if a == b:
                return 0
            try:
                if int(a) < int(b):
                    return -1
                else:
                    return 1
            except ValueError:
                return cmp(a, b)
        loci_list.sort(key=itemgetter(1))
        loci_list.sort(cmp=chr_cmp, key=itemgetter(0))

    return loci_set, loci_list


def parse_snp_chip_txt(filename):

    data = {}
    lineno = 0

    if filename.endswith(".gz"):
        from gzip import open as xopen
    else:
        xopen = open

    with xopen(filename) as f:
        lineno += 1
        hline = next(f).rstrip("\r\n")
        headers = hline.split("\t")[1:]
        for line in f:
            line = line.rstrip("\r\n")
            fields = line.split("\t")
            snpid = int(fields[0])
            data[snpid] = fields[1:]

            if len(data[snpid]) != len(headers):
                raise ParseError("Error: {}: line {}: mismatch between number of headers ({}) and data array length ({}).".format(filename, lineno, len(headers), len(data[snpid])))

    return headers, data


# SAM flag field
F_PAIRED        = 0x001 # the read is paired in sequencing
F_PAIR_MAPPED   = 0x002 # the read is mapped in a proper pair
F_UNMAPPED      = 0x004 # the query sequence itself is unmapped
F_MATE_UNMAPPED = 0x008 # the mate is unmapped
F_STRAND        = 0x010 # strand of the query (1 for reverse)
F_MATE_STRAND   = 0x020 # strand of the mate
F_FIRST_READ    = 0x040 # the read is the first read in a pair
F_SECOND_READ   = 0x080 # the read is the second read in a pair
F_SUBSEQUENT    = 0x100 # the alignment is not primary
F_QCFAIL        = 0x200 # QC failure
F_DUP           = 0x400 # optical or PCR duplicate
F_SUPP          = 0x800 # supplementary alignment

def parse_sam(filename):
    lineno = 0
    with open(filename) as f:
        for line in f:
            lineno += 1
            if line[0] == '@':
                # ignore headers
                continue
            line = line.rstrip("\r\n")
            fields = line.split("\t")

            if len(fields) < 11:
                raise ParseError("{}: line {}: expected at least 11 columns, got {}.".format(filename, lineno, len(fields)))

            qname = fields[0]   # query (pair) name
            flag = int(fields[1]) # bitwise flag (see above)
            rname = fields[2]   # reference sequence name
            pos = fields[3]     # 1-based lefmost position of clipped sequence
            mapq = fields[4]    # mapping quality (phred-scaled)
            cigar = fields[5]   # extended CIGAR string
            mrnm = fields[6]    # mate reference sequence name
                                # (`=' if same as rname)
            mpos = fields[7]    # 1-based mate position
            isize = fields[8]   # inferred insert size
            seq = fields[9]     # query sequence on the same strand as the ref
            qual = fields[10]   # query base quality scores (ASCII-33 phred)
            opt = fields[11:]   # optional fields in the format TAG:VTYPE:VAL
                                # (see mapper documentation for definitions)

            yield qname, flag, rname, pos, mapq, cigar, mrnm, mpos, isize, seq, qual, opt
