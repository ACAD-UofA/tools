#!/bin/sh

odir=stage3
odir1=${odir}/hetSNPs
odir2=${odir}/indels-etc

trap 'kill 0' EXIT

die() {
	echo Error: $@
	exit 1
}

do_filter() {
	vcf=$1
	sample=$2
	chr=$3

	# extract heterozygous SNPs
	bcftools view \
		-g het \
		-v snps \
		-O z \
		-o ${odir1}/${sample}.hetSNPs.${chr}.vcf.gz \
		${vcf} \
	|| die "${sample}:${chr}: bcftools view1"

	bcftools index ${odir1}/${sample}.hetSNPs.${chr}.vcf.gz \
	|| die "${sample}:${chr}: bcftools index1"

	# extract non-SNPS (indels, mnps, other)
	bcftools view \
		-V snps \
		-O z \
		-o ${odir2}/${sample}.indels-etc.${chr}.vcf.gz \
		${vcf} \
	|| die "${sample}:${chr}: bcftools view2"

	bcftools index ${odir2}/${sample}.indels-etc.${chr}.vcf.gz \
	|| die "${sample}:${chr}: bcftools index2"
}

export odir1
export odir2
export -f do_filter
export -f die

mkdir -p $odir1 $odir2

for vcf in stage2/*/tp-baseline.vcf.gz; do
	path=${vcf%/tp-baseline.vcf.gz}
	samplechr=$(basename $path)
	chr=${samplechr//*\./}
	sample=${samplechr%\.${chr}}
	echo do_filter $vcf $sample $chr
done | parallel -j 6
