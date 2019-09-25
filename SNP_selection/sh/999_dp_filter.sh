#!/bin/sh

numchrom=26 # sheep autosomes
odir=stage5/dp

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

do_dpstuff() {
	sample=$1

	python vcf_dp.py stage4/${sample}.raw_snps.*.vcf.gz \
		> ${odir}/${sample}.dp.txt
}

export odir
export -f die
export -f do_dpstuff

mkdir -p $odir

for s in $samples; do
	echo do_dpstuff $s
done | parallel -j 6
