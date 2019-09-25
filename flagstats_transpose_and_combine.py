#!/usr/bin/env python
"""
Collate the output from multiple runs of ``samtools flagstat foo.bam''.

Graham Gower, 2015.
"""

from __future__ import print_function
import sys
import re
import csv
import math
import os.path

if len(sys.argv) < 2:
    print("usage: {} file1.flagstats [file2.flagstats...]".format(sys.argv[0]))
    sys.exit(1)

re_line = re.compile(r"^([0-9]+) \+ ([0-9]+) (.*?)(?:\((.*?)%:(.*?)%\))?$")

# Extract column headers from the first file.
# Assumes all files are the same.
header = ["Sample",]
with open(sys.argv[1]) as f:
    lineno = 0
    for line in f:
        lineno += 1

        m = re_line.match(line)
        if m is None:
            print("{}: line {}: match failed.".format(sys.argv[1], lineno), file=sys.stderr)
            sys.exit(1)

        #n1 = m.group(1)
        #n2 = m.group(2)
        desc = m.group(3)

        if desc.startswith("in total"):
            desc = "total"

        desc_pass = "{} (QC pass)".format(desc)
        desc_fail = "{} (QC fail)".format(desc)

        header.append(desc_pass)
        header.append(desc_fail)

writer = csv.DictWriter(sys.stdout, fieldnames=header)
writer.writeheader()

for filename in sys.argv[1:]:
    with open(filename) as f:
        lineno = 0

        sample = os.path.basename(filename)[:-len(".flagstats")]
        row = {"Sample": sample}

        for line in f:
            lineno += 1

            m = re_line.match(line)
            if m is None:
                print("{}: line {}: match failed.".format(filename, lineno), file=sys.stderr)
                sys.exit(1)

            mg = m.groups()

            n1 = mg[0]
            n2 = mg[1]
            desc = mg[2]

            if mg[3]:
                try:
                    pc1 = float(mg[3])
                    if math.isnan(pc1):
                        raise ValueError()
                except ValueError:
                    pass
                else:
                    n1 += " ({}%)".format(pc1)

            if mg[4]:
                try:
                    pc2 = float(mg[4])
                    if math.isnan(pc2):
                        raise ValueError()
                except ValueError:
                    pass
                else:
                    n2 += " ({}%)".format(pc2)


            if desc.startswith("in total"):
                desc = "total"

            desc_pass = "{} (QC pass)".format(desc)
            desc_fail = "{} (QC fail)".format(desc)

            row[desc_pass] = n1
            row[desc_fail] = n2

    writer.writerow(row)

