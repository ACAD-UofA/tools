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
Covert snp array .txt to eigenstrat .geno/.ind/.snp files.
SNPs are sorted by SNPid.

e.g.
$ ./snpchip2eig.py \
        --bed /localscratch/BisonProjects/SNPs/9908_Bovid_SNP_10k.bed \
        -o 9908_Bovid749 \
        ~/ownCloud/Shared/Bison_projects/Bioinfo/40k_Data/Decoder/40843_749Bovid_Transposed_RefAltCode.txt.gz
"""

import sys
import random
import parsers

def eig_chr_map(chrom, nchrom=29):
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

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Convert SNPs from SNP chip .txt (with 0=REF/REF, 1=REF/ALT, 2=ALT/ALT, ?=unknown) to eigenstrat format (with 2=REF/REF, 1=REF/ALT, 0=ALT/ALT, 9=unknown).")
    parser.add_argument("--bed", required=True, help="SNP bed positions, containing SNP ids as fourth column")
    parser.add_argument("-o", "--oprefix", required=True, help="Output prefix for eigenstrat *.{geno,ind,snp} files")
    parser.add_argument("txt", help="SNP chip data (optionally gzipped)")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    headers, data = parsers.parse_snp_chip_txt(args.txt)

    with open("{}.ind".format(args.oprefix), "w") as f:
        for h in headers:
            hfields = h.split("_")

            # Come up with (semi) sensible group/sample names.
            if len(hfields) == 1:
                group = h
            else:
                if h[0] != "X":
                    group = "_".join(hfields[:-1])
                else:
                    if hfields[1][0] in ("A", "B"):
                        group = "_".join(hfields[2:])
                    else:
                        group = h
            if group.startswith("rep_"):
                group = group[4:]
            elif group.startswith("replig_"):
                group = group[7:]
            elif group.startswith("pool_"):
                group = group[5:]
            elif group.startswith("lig_col_"):
                group = group[8:]
            elif group.startswith("lig_eth_"):
                group = group[8:]
            elif group.startswith("3ul_"):
                group = group[4:]
            elif group.startswith("6ul_"):
                group = group[4:]
            elif group.startswith("10ul_"):
                group = group[5:]

            f.write("{}\tU\t{}\n".format(h, group))

    gt_map = {"0":"2", "1":"1", "2":"0", "?":"9"}
    #gt_map = {"0":"0", "1":"1", "2":"2", "?":"9"}

    _, loci_list = parsers.parse_bed(args.bed, sort_by_id=False)
    
    with open("{}.snp".format(args.oprefix), "w") as f_snp, \
         open("{}.geno".format(args.oprefix), "w") as f_geno:
        for chrom, pos, snpid in loci_list:
            f_snp.write("{}_{}\t{}\t0\t{}\n".format(chrom, pos, eig_chr_map(chrom), pos))

            gts = "".join(gt_map[gt] for gt in data[snpid])
            f_geno.write(gts)
            f_geno.write("\n")
