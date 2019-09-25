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

import sys
import parsers
import collections
import math

if __name__ == "__main__":

    counts = collections.defaultdict(lambda: collections.Counter())
    reads = set()
    lset = set()

    for line in parsers.parse_sam("/dev/stdin"):
        qname, flag, _, _, mapq, _, _, _, _, _, _, _ = line

        # filter reads
        if flag & (parsers.F_UNMAPPED|parsers.F_SUBSEQUENT|parsers.F_QCFAIL|parsers.F_SUPP) != 0:
            # ignore all these types of mappings
            continue

        # filter by mapping quality too
        try:
            if int(mapq) < 25:
               continue
        except ValueError:
            # Integer conversion failed. No mapping quality provided.
            continue

        q_fields = qname.split(":")
        read_id = ":".join(q_fields[:7])

        if len(q_fields) == 7:
            # master read
            if read_id in reads:
                raise Exception("Error: multiple mappings for master read {}".format(qname))
            reads.add(read_id)
            continue

        rep = int(q_fields[-2])
        length = int(q_fields[-1])

        counts[read_id][length] += 1
        lset.add(length)

    lengths = list(lset)
    lengths.sort(reverse=True)

    print("#RID\t" + "\t".join(str(l) for l in lengths))

    for read_id, counter in counts.iteritems():
        counts = "\t".join(str(counter[l]) for l in lengths)
        print("{}\t{}".format(read_id, counts))
