#!/bin/sh

vcfdir=/localscratch/grg/bovid_genomes/OARv40/vcf/
suffix=.FILTER_LowQual.vcf.gz
odir=het_with_ti
autosomes=26

vcfhet1() {
	vcf=$1
	sample=$2

	./vcfhet -t -r X $vcf > ${odir}/${sample}.X
	./vcfhet -t -r MT $vcf > ${odir}/${sample}.MT

	autlist=$(seq -s, $autosomes)
	./vcfhet -t -r $autlist $vcf > ${odir}/${sample}.Aut
}

vcfhet2() {
	vcf=$1
	sample=$2

	./vcfhet -r X $vcf > ${odir}/${sample}.X
	./vcfhet -r MT $vcf > ${odir}/${sample}.MT

	autlist=$(seq -s, $autosomes)
	./vcfhet -r $autlist $vcf > ${odir}/${sample}.Aut
}

export -f vcfhet2
export odir
export autosomes

mkdir -p $odir

for vcf in $vcfdir/*$suffix; do
	sample=$(basename ${vcf%$suffix})
	echo vcfhet2 $vcf $sample
done | parallel -j 8

