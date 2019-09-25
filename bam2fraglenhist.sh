#!/bin/sh

module load zlib samtools

# 0xf01 == 0x1|0x100|0x200|0x400|0x800
samtools view \
	$1 \
	| awk ' and($2, 0xf01)==0 { l[length($10)]++ }
		END {
			for (i=0; i<500; i++)
				printf("%d\t%d\n", i, l[i])
		}'
