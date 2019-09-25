#!/usr/bin/env python
"""
Accept a .geno file as input and output a new .geno file which specifies
1 for missing data, 0 for data present.
"""

import sys

for line in sys.stdin:
    line = line.rstrip("\r\n")
    new_line = []
    for gt in line:
        if gt == '9':
            new_line.append("1")
        else:
            new_line.append("0")
    sys.stdout.write("".join(new_line))
    sys.stdout.write("\n")
