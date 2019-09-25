#!/bin/sh

ref=/localscratch/Refs/Ovis_aries/Oar_v4.0/OARv40.fasta
intsz1=60 # interval size, each side of the site of interest
numchrom=26 # sheep autosomes
odir=stage6

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

do_winlist() {
	vcf=$1
	chr=$2
	intsz=$3
	chrlen=$(awk '$1=='${chr}' {print $2}' ${ref}.fai)

	zcat $vcf \
	 | awk '$1=='${chr}' {
			from=$2-'${intsz}'
			if (from <= 0) from = 1
			to=$2+'${intsz}'
			if (to > '${chrlen}') to = '${chrlen}'
			printf "%s:%d-%d\n", $1, from, to
		}' \
	 | xargs samtools faidx ${ref}
}

do_cnswin() {
	sample=$1
	chr=$2
	cnsvcf=stage2/${sample}.${chr}/tp-baseline.vcf.gz
	hetvcf=stage5/${sample}.snps-no-indels.${chr}.vcf.gz
	opfx1=${odir}/${sample}.snpwin${intsz1}.${chr}

	# generate a list of chr:from-to windows from the het vcf
	do_winlist ${hetvcf} ${chr} ${intsz1} \
		| bcftools consensus ${cnsvcf} \
		> ${opfx1}.fa 2>/dev/null \
	|| die "${sample}:${chr} bcftools consensus"
}

export ref
export odir
export intsz1
export -f die
export -f do_winlist
export -f do_cnswin

mkdir -p ${odir}

#echo do_cnswin A16171_Bpriscus_Yukon_40999 25 | parallel

samples="
Ovis_aries_0
"

for s in $samples; do
	for n in `seq $numchrom`; do
		echo do_cnswin $s $n
	done
done | parallel -j 8

