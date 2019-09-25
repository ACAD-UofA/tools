#!/usr/bin/env python
# Copyright (c) 2016 Graham Gower <graham.gower@gmail.com>
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
Convert C->T in a fastq file, ignoring cytosines in a CpG context.
Outputs converted forward sequence and converted reverse sequence.
"""

from __future__ import print_function
import sys
import parsers

def revcomp(seq):
    """
    Reverse complement of sequence @seq.
    """
    revmap = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
    return "".join((revmap[s] for s in reversed(seq)))

def bstreat(seq):
    """
    Convert C->T, ignoring cytosines in a CpG context.
    """
    s = []
    for i in xrange(len(seq)):
        if seq[i] != "C":
            s.append(seq[i])
        else:
            if i==0:
                # start of the sequence
                if seq[i+1] == "G":
                    s.append("C")
                else:
                    s.append("T")
                continue
            elif i==len(seq)-1:
                # end of the sequence
                if seq[i-1] == "G":
                    s.append("C")
                else:
                    s.append("T")
                continue
            else:
                if seq[i+1] == "G" or seq[i-1] == "G":
                    s.append("C")
                else:
                    s.append("T")
                continue

    return "".join(s)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: {} file.fastq[.gz]".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    filename = sys.argv[1]
    
    for lineno, line in enumerate(parsers.parse_fa(filename),1):

        if len(line) == 3:
            # fastq, output single sequence
            label, seq, qual = line
            bs_seq = bstreat(seq)
            print("{}\n{}\n+\n{}".format(label, bs_seq, qual))

        elif len(line) == 2:
            # fasta, output two sequences
            label, seq = line
            lfields = label.split()
            sname = lfields[0]
            if len(lfields) > 1:
                comment = " " + " ".join(lfields[1:])
            else:
                comment = ""

            # strip trailing pipe, if it exists
            if sname[-1] == "|":
                sname = sname[:-1]

            bs_seq = bstreat(seq)
            print("{}|BS|F{}\n{}\n".format(sname, comment, bs_seq))
            bs_seq = None # allow garbage collector to reclaim memory
            bs_rseq = bstreat(revcomp(seq))
            print("{}|BS|R{}\n{}\n".format(sname, comment, bs_rseq))
