#!/bin/sh
# Filter out het sites with proximity to an indel or other non-SNP variant

intsz=60 # interval size, each side of the site of interest
numchrom=26 # sheep autosomes
odir=stage5

samples="
Ovis_aries_0
"

trap 'kill 0' EXIT

die() {
	echo Error: $@
	exit 1
}

do_filter() {
	sample=$1
	chr=$2
	hetvcf=stage3/hetSNPs/${sample}.hetSNPs.${chr}.vcf.gz
	excvcf=stage3/indels-etc/${sample}.indels-etc.${chr}.vcf.gz
	ovcf=${odir}/${sample}.snps-no-indels.${chr}.vcf.gz

	# generate a list of windows for exclusion
	excwinlist=stage3/indels-etc/${sample}.indels-etc.${chr}.bed
	zcat ${excvcf} \
	 | awk '$1=='${chr}' {
			l4 = length($4)
			l5 = length($5)
			if (l4>l5)
				varlen = l4
			else
				varlen = l5
			from=$2 - varlen - '${intsz}' -1
			to=$2 + varlen + '${intsz}'
			printf "%s\t%d\t%d\n", $1, from, to
		}' \
	 > ${excwinlist}

	bcftools view -T ^${excwinlist} \
		-O z \
		-o ${ovcf} \
		${hetvcf} \
	|| die "${sample}:${chr}: bcftools view"

	bcftools index ${ovcf} \
	|| die "${sample}:${chr}: bcftools index"
}

export odir
export intsz
export -f die
export -f do_filter

mkdir -p ${odir}

for s in $samples; do
	for n in `seq $numchrom`; do
		echo do_filter $s $n
	done
done | parallel -j 8
