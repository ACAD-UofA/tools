#!/bin/sh

samples="
3A_BisonXUrals
5A_BisonXUrals
6A_BisonXUrals
7A_BisonXUrals
17A_BisonXUrals
18A_BisonXUrals
4093A_Bisonbonasus
"

sampledir="/localscratch/jsoubrier/150623_IMS_BLl_BHC04Bison10K/2_Paleomix/CowRef"
vcfoutdir=./
ref="/localscratch/Refs/Bos_taurus/Bos_taurus_NoUnplaced_BWA6_2_2013_04/REFBTaurusNoUnplaced.fasta"
bed="/localscratch/BisonProjects/SNPs/9908_Bovid_SNP_10k.bed"

module load zlib/1.2.8-gnu_4.8.0
module load samtools/1.2
module load bcftools/1.2

for s in $samples; do
	bam_in="${sampledir}/${s}.REFBTaurusNoUnplaced.realigned.bam"
	vcf_out="${vcfoutdir}/${s}_10k_snps.vcf.gz"

	if [ -f "$vcf_out" ] && [ x"$1" != x"-f" ]; then
		continue
	fi

	#
	# 1) Pile up aligned reads, calculate likelihoods (the bottleneck)
	# 2) Call genotypes
	#
	samtools mpileup \
		-f $ref \
		-l $bed \
		-t DP,DPR \
		--skip-indels \
		-u \
		$bam_in \
	| \
	bcftools call \
		-c \
		-O z \
		-o $vcf_out \
	&
done

wait

