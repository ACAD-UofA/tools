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
Print the nucleotide composition of the specified fasta or fastq file.
"""

from __future__ import print_function
import sys
import parsers
import collections

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: {} in.{{fa,fq}}".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    ntcomp = collections.Counter()

    for seq in parsers.parse_fa(sys.argv[1]):
        s = seq[1]
        ntcomp.update(s)

    gc = float(ntcomp['G']+ntcomp['C']) / (ntcomp['G']+ntcomp['C']+ntcomp['A']+ntcomp['T'])
    print(ntcomp)
    print("GC ratio = {:.3f}".format(gc))
