#!/bin/sh

idir=/localscratch/grg/10k/BisonX/02_all_10k_eig/
pfx=`cat /localscratch/grg/10k/BisonX/prefix`
f_geno=${idir}/${pfx}.geno
f_ind=${idir}/${pfx}.ind

awk '
BEGIN {OFS="\t"}
{
	if ($3 == "AncientWisentA" || $3 == "HistoricalWisentA" || $3 == "ModernWisentA")
		$3 = "WisentA"
	print
}' ${f_ind} > WisentA_vs_WisentB.ind

./f-stats.py \
	--geno $f_geno \
	--ind $f_ind \
	--groups \
	WisentB,AncientWisentA,Steppe,Bubalus_bubalis \
	WisentB,HistoricalWisentA,Steppe,Bubalus_bubalis \
	WisentB,ModernWisentA,Steppe,Bubalus_bubalis \
	> branches.tsv

./f-stats.py \
	--geno $f_geno \
	--ind WisentA_vs_WisentB.ind \
	--groups \
	WisentB,WisentA,Steppe,Bubalus_bubalis \
	| tail -n +2 >> branches.tsv

exit

./f-stats.py \
	--geno $f_geno \
	--ind $f_ind \
	--groups \
	WisentB_6A,AncientWisentA_4093A,Steppe_A875,Bubalus_bubalis \
	WisentB_6A,HistoricalWisentA_13947A,Steppe_A875,Bubalus_bubalis \
	WisentB_6A,ModernWisentA,Steppe_A875,Bubalus_bubalis \
	WisentB_6A,AncientWisentA_4093A,Steppe_A3133,Bubalus_bubalis \
	WisentB_6A,HistoricalWisentA_13947A,Steppe_A3133,Bubalus_bubalis \
	WisentB_6A,ModernWisentA,Steppe_A3133,Bubalus_bubalis \
	> ry_branches.tsv

