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
Convert snp array .txt to phylip format. SNPs are sorted by SNPid

e.g.
$ ./snpchip2phylip.py \
        --bed /localscratch/BisonProjects/SNPs/9908_Bovid_SNP_10k.bed \
        ~/ownCloud/Shared/Bison_projects/Bioinfo/40k_Data/Decoder/40843_749Bovid_Transposed_RefAltCode.txt.gz \
        > 9908_749Bovid.phy
"""

import sys
import parsers

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Convert SNPs from SNP chip .txt to phylip format. (both with 0=REF/REF, 1=REF/ALT, 2=ALT/ALT, ?=unknown)")
    parser.add_argument("--bed", required=True, help="SNP bed positions, containing SNP ids as fourth column")
    parser.add_argument("txt", help="SNP chip data (optionally gzipped)")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    headers, data = parsers.parse_snp_chip_txt(args.txt)
    _, loci_list = parsers.parse_bed(args.bed)
 
    sys.stdout.write(" {} {}\n".format(len(headers), len(loci_list)))

    for i, h in enumerate(headers):
        sys.stdout.write("{} ".format(h))
        gts = "".join(data[snpid][i] for _ ,_, snpid in loci_list)
        sys.stdout.write(gts)
        sys.stdout.write("\n")
