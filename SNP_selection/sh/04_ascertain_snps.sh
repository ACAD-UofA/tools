#!/bin/sh

numchrom=26 # sheep autosomes
odir=stage4

#TODO: A15654_Bcaucasicus_CaucasusKubanOblast_39

samples="
A16121_Bpriscus_Yukon_14792
A16171_Bpriscus_Yukon_40999
A1_BisonX_UralsRasik1_14874
A2494_Bprimigenius_RussiaYenseiRiver_14082
A3020_Bpriscus_Minnesota_7642
A3133_Bpriscus_Yukon_30612
A860_Bpriscus_Alaska_33154
A875_Bpriscus_Siberia_50k
BBO3569_Bbonasus_Bialowieza_0
BBO3574_Bbonasus_Bialowieza_0
BB_20087_Bbison_bison_0
Q229_Bbonasus_0
SW18_Bbison_athabascae_0
"

trap 'kill 0' EXIT

die() {
	echo Error: $@
	exit 1
}

do_ascertain() {
	sample=$1
	chr=$2

	bcftools isec \
		-r ${chr} \
		-c all \
		-i 'GT=="0/0"' \
		-i 'GT=="0/1"' \
		-O z \
		-o ${odir}/${sample}.raw_snps.${chr}.vcf.gz \
		-n=2 \
		-w2 \
		stage3/hom/hom_consensus.${chr}.vcf.gz \
		stage3/hetSNPs/${sample}.hetSNPs.${chr}.vcf.gz \
	|| die "${sample}:${chr}: bcftools isec"

	bcftools index ${odir}/${sample}.raw_snps.${chr}.vcf.gz \
	|| die "${sample}:${chr}: bcftools index"
}

export odir
export -f die
export -f do_ascertain

mkdir -p $odir

for s in $samples; do
	for n in `seq $numchrom`; do
		echo do_ascertain $s $n
	done
done | parallel -j 6
