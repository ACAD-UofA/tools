#!/bin/sh

numchrom=26 # sheep autosomes
odir=02_depth

trap 'kill 0' EXIT

die() {
	echo Error: $@
	exit 1
}

do_filter() {
	sample=$1
	chr=$2
	shift; shift
	dp=$@

	bed=${odir}/${sample}.snpwin60.${chr}.bed
	ivcf=01_WindowMasker/${sample}.snpwin60.${chr}.vcf.gz
	ovcf=${odir}/${sample}.snpwin60.${chr}.vcf.gz

	python filter_depth.py $ivcf $dp > $bed

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

mkdir -p $odir

samples="
Ovis_aries_0
"

dp_dir=$HOME/src/grg-utils/src/vcfhet/het/
thres=0.5 # top 50%
for sample in $samples; do
	dplist=$(awk '{print $1; if ($3>'$thres') exit}' ${dp_dir}/${sample}.DP)
	for chr in `seq $numchrom`; do
		echo do_filter $sample $chr $dplist
	done
done | parallel -j 8
