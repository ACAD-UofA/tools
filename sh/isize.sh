#!/bin/sh

bam=$1
nchr=26

#echo -n "sample"
#for c in `seq $nchr` X; do
#	echo -e "\t$c"
#done

echo -n $(basename ${bam%.bam})

for c in `seq $nchr` X; do
	samtools view \
		-q 25 \
		-F `python -c 'print(hex(0x4|0x100|0x200|0x400|0x800))'` \
		$bam \
		chr$c \
	| awk '{
		n++;
		l[n] = length($10);
		sum += l[n];
	}
	END {
		if (n > 1) {
			mean = sum/n;
			for (i=1; i<=n; i++)
				sumsq += (mean-l[n])*(mean-l[n]);
			var = sumsq/(n-1);
		} else {
			mean = var = 0;
		}
		printf "\t%d\t%.3lf\t%.3lf", n, mean, sqrt(var);
	}
	'
done

echo
