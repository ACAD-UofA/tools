#!/bin/sh

grep -A1 ^@ $1 \
	| awk '!/^@/ && !/^--/ { l[length($1)]++ }
		END {
			for (i=0; i<500; i++)
				printf("%d\t%d\n", i, l[i])
		}'
