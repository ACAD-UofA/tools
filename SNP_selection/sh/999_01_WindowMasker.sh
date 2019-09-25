#!/bin/sh

numchrom=26 # sheep autosomes
odir=01_WindowMasker

trap 'kill 0' EXIT

die() {
	echo Error: $@
	exit 1
}

do_filter() {
	sample=$1
	chr=$2

	fa_in=../stage6/${sample}.snpwin60.${chr}.fa
	fa_out=${odir}/${sample}.snpwin60.${chr}.fa
	bed=${odir}/${sample}.snpwin60.${chr}.bed
	ivcf=../stage5/${sample}.snps-no-indels.${chr}.vcf.gz
	ovcf=${odir}/${sample}.snpwin60.${chr}.vcf.gz

	python filter_WindowMasker.py $fa_in > $fa_out

	for pos in $(grep "^>" $fa_out | cut -f2 -d: | cut -f1 -d-); do
		from=$((pos+60-1))
		to=$((pos+60))
		echo -e "$chr\t$from\t$to"
	done > $bed

	bcftools view -T $bed \
		-O z \
		-o ${ovcf} \
		${ivcf} \
	|| die "${sample}:${chr}: bcftools view"

	bcftools index ${ovcf} \
	|| die "${sample}:${chr}: bcftools index"
}

export odir
export -f die
export -f do_filter

samples="
Ovis_aries_0
"

for sample in $samples; do
	for chr in `seq $numchrom`; do
		echo do_filter $sample $chr
	done
done | parallel -j 8
