#!usr/bin/env python

import pysam

with pysam.Samfile("-", "rb") as infile:
    with pysam.Samfile("-", "wb") as outfile:
        for (read_num, read) in enumerate(infile):
            if read.flag in (0x63, 0x93, 0x53, 0xa3):
                outfile.write(read)
