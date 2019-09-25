#!/bin/sh
# Determine sex using number of reads mapping to X and autosome as a
# proxy for gene dosage.  This follows Skoglund et al. [1], but uses
# X and autosome instead of X and Y.
# [1] Skoglund et al. 2013. doi:10.1016/j.jas.2013.07.004
 

get_sex() {
	bam=$1

	if [ ! -r ${bam}.bai -a ! -r ${bam%.bam}.bai ]; then
		echo "Error: could not find ${bam}.bai, nor ${bam%.bam}.bai. Try running:"
		echo " samtools index ${bam}"
		exit 1
	fi

	file=$(basename ${bam%.bam})

	samtools idxstats $bam \
	| awk '
		$1~/^([Cc]hr)?[1-9][0-9]?$/ {
			La += $2
			Na += $3
		}
		$1~/^([Cc]hr)?X$/ {
			Lx = $2
			Nx = $3
		}
		END {
			if (Nx+Na == 0) { exit 1 }
			Rx = Nx/(Nx+Na)
			ERx = Lx/(Lx+La) # Expectation from the lengths
			SE = sqrt(Rx*(1-Rx)/(Nx+Na))
			printf "'$file'\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf-%lf\n",
				Nx, Lx, Na, La, ERx, 0.5*ERx, Rx, SE, Rx-1.96*SE, Rx+1.96*SE
		}'
}

if [ $# -lt 2 ]; then
	echo "usage: $0 file1.bam [... fileN.bam]"
	exit 1
fi

echo -e "File\tNx\tLx\tNa\tLa\tERx\t0.5*ERx\tRx\tSE[Rx]\t95%CI[Rx]"
for b in "$@"; do
	get_sex $b
done
