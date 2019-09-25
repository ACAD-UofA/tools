#!/usr/bin/env python
# Randomize the order of the fasta sequences given on stdin.  Sequences are
# given new names based on their order in the output.
#
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

import sys
import random
from parsers import parse_fa

sequences = []

for line in parse_fa("/dev/stdin"):
    assert(len(line) == 2)
    _, seq = line
    sequences.append(seq)

random.shuffle(sequences)

def print_fa(label, seq, linelen=60):
    print(">{}".format(label))
    while len(seq) > 0:
        print(seq[:linelen])
        seq = seq[linelen:]

for n, seq in enumerate(sequences, 1):
    print_fa(n, seq)
