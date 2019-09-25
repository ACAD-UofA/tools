#!/bin/sh

numchrom=26 # sheep autosomes
odir=stage3/hom

trap 'kill 0' EXIT

die() {
	echo Error: $@
	exit 1
}

do_hom_consensus() {
	n=$1

	bcftools isec \
		-r ${n} \
		-c all \
		-i 'GT=="0/0"' \
		-O z \
		-o ${odir}/hom_consensus.${n}.vcf.gz \
		-n=8 \
		-w1 \
		vcf_samtools_bcftools/Ovis_aries_0.FILTER_LowQual.vcf.gz \
		vcf_samtools_bcftools/Giraffa_camelopardalis_0.FILTER_LowQual.vcf.gz \
		vcf_unified_genotyper/Ovis_aries_0.${n}.vcf.gz \
		vcf_unified_genotyper/Gcamelopardalis_KenyaMA1_0.${n}.vcf.gz \
		vcf_unified_genotyper/Gcamelopardalis_NashvilleZOO_0.${n}.vcf.gz \
		vcf_haplotype_caller/Ovis_aries_0.${n}.vcf.gz \
		vcf_haplotype_caller/Gcamelopardalis_KenyaMA1_0.${n}.vcf.gz \
		vcf_haplotype_caller/Gcamelopardalis_NashvilleZOO_0.${n}.vcf.gz \
	|| die "${sample}:${n}: bcftools isec"

	bcftools index ${odir}/hom_consensus.${n}.vcf.gz \
		|| die "${sample}:${n}: bcftools index"
}

export odir
export -f die
export -f do_hom_consensus

mkdir -p $odir

for n in `seq $numchrom`; do
	echo do_hom_consensus $n
done | parallel -j 6
