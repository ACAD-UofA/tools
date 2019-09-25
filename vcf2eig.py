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
Produce SNP data from one or more VCFs for use with Eigensoft.
SNPs are sorted by SNPid.

e.g.
$ ./vcf2eig.py \
        --min-depth 2 \
        --ind-group-label ANCIENT \
        --bed /localscratch/BisonProjects/SNPs/9908_Bovid_SNP_10k.bed \
        -o 9908_Ancient_Bison \
        /localscratch/jsoubrier/10k/BisonX/run_raxml/*.vcf.gz
"""

from __future__ import print_function
import sys
import collections
import random
import parsers

def eig_chr_map(chrom, nchrom=29):
    """
    Translate chromosome names, according to Eigensoft's weird wants.
    """
    if chrom == "X":
        return nchrom+1
    elif chrom == "Y":
        return nchrom+2
    elif chrom in ("M","MT","Mt"):
        return 90
    else:
        try:
            chrom_int = int(chrom)
        except ValueError:
            raise RuntimeError("Unknown chromosome {}".format(chrom_str))
        return chrom_int


def create_eig(loci_list, variants, samples, prefix, exclude_missing, snp_6col, ind_group_label, ignore_het, haploidisation, ignore_transitions):
    """
    Produce output files in the EIGENSTRAT format expected by Eigensoft.

    Format specifications taken from Eigensoft 6.0.1: CONVERTF/README.
    """

    genotype_file = "{}.geno".format(prefix)
    snp_file = "{}.snp".format(prefix)
    ind_file = "{}.ind".format(prefix)

    # track snp_ids for which there is data for all samples
    idlist_complete = set()

    """
    The genotype file contains 1 line per SNP.  
      Each line contains 1 character per individual:
      0 means zero copies of reference allele.
      1 means one copy of reference allele.
      2 means two copies of reference allele.
      9 means missing data.
    """
    with open(genotype_file, "w") as f:
        for chrom, pos, snp_id in loci_list:
            gt_list = []
            exclude = False
            for sample in samples:
                try:
                    _, _, ref, alt, gt_str = variants[sample][(chrom, pos)]
                except KeyError:
                    if exclude_missing:
                        exclude = True
                        break
                    gt = 9
                else:
                    # GT field:
                    #   0/0 = REF/REF
                    #   0/1 = REF/ALT1
                    #   1/1 = ALT1/ALT1
                    #   1/2 = ALT1/ALT2
                    gt_fields = gt_str.split("/")

                    if "." in gt_fields:
                        gt_list.append(9)
                        continue

                    if ignore_transitions:
                        alt_i = int(gt_fields[1]) -1
                        if alt_i >= 0:
                            alt_fields = alt.split("/")
                            alt = alt_fields[alt_i]
                            if ((ref == "C" and alt == "T") or
                               (ref == "T" and alt == "C") or
                               (ref == "A" and alt == "G") or
                               (ref == "G" and alt ==  "A")):
                                gt_list.append(9)
                                continue

                    if ignore_het and gt_fields[0] != gt_fields[1]:
                        gt = 9
                    elif haploidisation and gt_fields[0] != gt_fields[1]:
                        gt = random.choice((0, 2))
                    else:
                        gt = 0
                        for x in gt_fields:
                            if x == "0":
                                gt += 1

                gt_list.append(gt)

            if exclude:
                continue

            idlist_complete.add(snp_id)

            for gt in gt_list:
                f.write(str(gt))
            f.write("\n")

    """
    The snp file contains 1 line per SNP.  There are 6 columns (last 2 optional):
      1st column is SNP name
      2nd column is chromosome.  X chromosome is encoded as 23.
        Also, Y is encoded as 24, mtDNA is encoded as 90, and XY is encoded as 91.
        Note: SNPs with illegal chromosome values, such as 0, will be removed
      3rd column is genetic position (in Morgans).  If unknown, ok to set to 0.0.
      4th column is physical position (in bases)
      Optional 5th and 6th columns are reference and variant alleles.
        For monomorphic SNPs, the variant allele can be encoded as X (unknown).
    """
    with open(snp_file, "w") as f:
        for chrom, pos, snp_id in loci_list:
            if exclude_missing and snp_id not in idlist_complete:
                continue

            if not snp_6col:
                # 4 column version.
                f.write("{}\t{}\t0.0\t{}\n".format(snp_id, eig_chr_map(chrom), pos))
                continue

            # 6 column version.
            #
            # Need to obtain consensus REF/ALT from all samples,
            # where REF might still be unknown.
            # This probably makes more sense with --exclude-missing
            #
            ref = None
            alt = None
            for sample in samples:
                if (chrom, pos) not in variants[sample]:
                    continue

                _, _, ref2, alt2, _ = variants[sample][(chrom, pos)]

                if ref and ref != ref2:
                    print("Warning: REF disagreement for SNP {} at locus {}:{}".format(snp_id, chrom, pos), file=sys.stderr)

                ref = ref2

                if alt2 != ".":
                    if alt and alt != alt2:
                        # Don't know who to believe...
                        alt = "X"
                    else:
                        alt = alt2

            if ref is None:
                # This is very noisy without --exclude-missing
                #print("Warning: REF for SNP {} at {}:{} is not known".format(snp_id, chrom, pos), file=sys.stderr)
                # XXX: What are the implications of this?
                ref = "X"

            if alt is None:
                alt = "X"

            f.write("{}\t{}\t0.0\t{}\t{}\t{}\n".format(snp_id, eig_chr_map(chrom), pos, ref, alt))

    """
    The indiv file contains 1 line per individual.  There are 3 columns:
      1st column is sample ID.  Length is limited to 39 characters, including
        the family name if that will be concatenated.
      2nd column is gender (M or F).  If unknown, ok to set to U for Unknown.
      3rd column is a label which might refer to Case or Control status, or
        might be a population group label.  If this entry is set to "Ignore", 
        then that individual and all genotype data from that individual will be
        removed from the data set in all convertf output.
    """
    with open(ind_file, "w") as f:
        for sample in samples:
            if ind_group_label:
                group = ind_group_label
            else:
                group = sample
            f.write("{}\tU\t{}\n".format(sample[:39], ind_group_label))

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Produce SNP data from one or more VCFs for use with Eigensoft (with 2=REF/REF, 1=REF/ALT, 0=ALT/ALT, 9=unknown)")
    parser.add_argument("--exclude-missing", action="store_true", default=False, help="Exclude SNPs lacking genotype data for any sample in the vcf [%(default)s]")
    parser.add_argument("--snp-6col", action="store_true", default=False, help="Produce a 6 column .snp file, with REF/ALT fields instead of 4 columns. [%(default)s]")
    parser.add_argument("--min-depth", type=int, default=5, help="Filter calls based upon the DP field [%(default)s]")
    parser.add_argument("--min-qual", type=int, default=0, help="Filter calls based upon the QUAL field [%(default)s]")
    parser.add_argument("--ind-group-label", default=None, help="Common group label for all samples for the .ind file")
    parser.add_argument("--ignore-het", action="store_true", default=False, help="Assume REF/ALT calls are unknown [%(default)s]")
    parser.add_argument("--haploidisation", action="store_true", default=False, help="Do not call hetereozygotes. Randomly select an allele.")
    parser.add_argument("--ignore-transitions", action="store_true", default=False, help="Set transitions (A<->G or C<->T) as missing data.")
    parser.add_argument("--bed", help="SNP bed positions, containing SNP ids as fourth column")
    parser.add_argument("-o", "--oprefix", required=True, help="Output prefix for eigenstrat *.{geno,ind,snp} files")
    parser.add_argument("vcf", nargs="+", help="One or more VCFs containing genotype data (optionally gzipped)")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if args.ignore_het and args.haploidisation:
        print("Error: --ignore-het and --haploidisation parameters are mutually exclusive", file=sys.stderr)
        exit(1)

    loci_set, loci_list = parsers.parse_bed(args.bed, sort_by_id=False, sort_by_coord=False)

    all_samples = []
    variants = collections.defaultdict(dict)

    for filename in args.vcf:
        vp = parsers.vcf_parser(filename)
        samples = next(vp)
        all_samples.extend(samples)
        for vcf_line in vp:
            (lineno, line, chrom, pos, id, ref, alt, qual, filter, info_str, fmt_fields, genotypes) = vcf_line

            if (chrom, pos) not in loci_set:
                continue

            if qual < args.min_qual:
                continue

            for sample in samples:
                gt = genotypes[sample]["GT"]
                if int(genotypes[sample]["DP"]) < args.min_depth:
                    gt = "./."

                variants[sample][(chrom, pos)] = (chrom, pos, ref, alt, genotypes[sample]["GT"])

    create_eig(loci_list, variants, all_samples, args.oprefix, args.exclude_missing, args.snp_6col, args.ind_group_label, args.ignore_het, args.haploidisation, args.ignore_transitions)
