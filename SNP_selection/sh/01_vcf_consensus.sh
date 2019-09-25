#!/bin/sh

sdf=/localscratch/Refs/Ovis_aries/Oar_v4.0/OARv40.sdf
numchrom=26 # sheep autosomes

samples="
A15654_Bcaucasicus_CaucasusKubanOblast_39
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
Gcamelopardalis_KenyaMA1_0
Gcamelopardalis_NashvilleZOO_0
Ovis_aries_0
Q229_Bbonasus_0
SW18_Bbison_athabascae_0
"

trap 'kill 0' EXIT

die() {
	echo Error: $@
	exit 1
}

do_consensus() {
	sample=$1
	n=$2

	vcf1=vcf_samtools_bcftools/${sample}.FILTER_LowQual.vcf.gz
	vcf2=vcf_unified_genotyper/${sample}.${n}.vcf.gz
	vcf3=vcf_haplotype_caller/${sample}.${n}.vcf.gz
	odir1=stage1/${sample}.${n}
	odir2=stage2/${sample}.${n}

	#mkdir -p ${odir1}

	export RTG_MEM="10g"
	export RTG_JAVA_OPTS="-XX:ParallelGCThreads=4"

	if [ ! -f ${odir1}/tp-baseline.vcf.gz.tbi ]; then
		rtg vcfeval \
			-b ${vcf1} \
			-c ${vcf2} \
			-t ${sdf} \
			--region=${n} \
			-o ${odir1} \
		|| die "vcfeval:1:${sample}:${n}"
	fi

	#mkdir -p ${odir2}

	rtg vcfeval \
		-b ${odir1}/tp-baseline.vcf.gz \
		-c ${vcf3} \
		-t ${sdf} \
		--region=${n} \
		-o ${odir2} \
	|| die "vcfeval:2:${sample}:${n}"
}

export sdf
export -f die
export -f do_consensus

mkdir -p stage1 stage2

for sample in $samples; do
	#for n in `seq $numchrom`; do
	for n in `seq 25`; do
		echo do_consensus $sample $n
	done
done | parallel -j 4
