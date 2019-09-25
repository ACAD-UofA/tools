#!/bin/sh

idir=/localscratch/grg/10k/BisonX/02_all_10k_eig/
pfx=`cat /localscratch/grg/10k/BisonX/prefix`
f_geno=${idir}/${pfx}.geno
f_ind=${idir}/${pfx}.ind
f_snp=${idir}/${pfx}.snp

awk '
BEGIN {OFS="\t"}
{
	if ($3 == "AncientWisent" || $3 == "HistoricalWisent" || $3 == "ModernWisent")
		$3 = "Wisent"
	print
}' ${f_ind} > Wisent_vs_CladeX.ind

./aabb.py \
	--geno $f_geno \
	--ind $f_ind \
	--snp $f_snp \
	--groups \
		AncientWisent,Steppe,Aurochs,Bubalus_bubalis \
		CladeX,Steppe,Aurochs,Bubalus_bubalis \
		HistoricalWisent,Steppe,Aurochs,Bubalus_bubalis \
		ModernWisent,Steppe,Aurochs,Bubalus_bubalis \
	> aabb-hypergeometric.tsv

./aabb.py \
	--geno $f_geno \
	--ind $f_ind \
	--snp $f_snp \
	--groups \
		AncientWisent,Aurochs,Steppe,Bubalus_bubalis \
		CladeX,Aurochs,Steppe,Bubalus_bubalis \
		HistoricalWisent,Aurochs,Steppe,Bubalus_bubalis \
		ModernWisent,Aurochs,Steppe,Bubalus_bubalis \
	| tail -n +2 >> aabb-hypergeometric.tsv

./aabb.py \
	--geno $f_geno \
	--ind Wisent_vs_CladeX.ind \
	--snp $f_snp \
	--groups \
		Wisent,Steppe,Aurochs,Bubalus_bubalis \
		CladeX,Steppe,Aurochs,Bubalus_bubalis \
	| tail -n +2 >> aabb-hypergeometric.tsv

./aabb.py \
	--geno $f_geno \
	--ind Wisent_vs_CladeX.ind \
	--snp $f_snp \
	--groups \
		Wisent,Aurochs,Steppe,Bubalus_bubalis \
		CladeX,Aurochs,Steppe,Bubalus_bubalis \
	| tail -n +2 >> aabb-hypergeometric.tsv

exit 1

#
# Consensus version

./aabb.py \
	--geno $f_geno \
	--ind $f_ind \
	--consensus \
	--snp $f_snp \
	--groups \
		AncientWisent,Steppe,Aurochs,Bubalus_bubalis \
		CladeX,Steppe,Aurochs,Bubalus_bubalis \
		HistoricalWisent,Steppe,Aurochs,Bubalus_bubalis \
		ModernWisent,Steppe,Aurochs,Bubalus_bubalis \
	> aabb-hypergeometric-cns.tsv

./aabb.py \
	--geno $f_geno \
	--ind $f_ind \
	--consensus \
	--snp $f_snp \
	--groups \
		AncientWisent,Aurochs,Steppe,Bubalus_bubalis \
		CladeX,Aurochs,Steppe,Bubalus_bubalis \
		HistoricalWisent,Aurochs,Steppe,Bubalus_bubalis \
		ModernWisent,Aurochs,Steppe,Bubalus_bubalis \
	| tail -n +2 >> aabb-hypergeometric-cns.tsv

./aabb.py \
	--geno $f_geno \
	--ind Wisent_vs_CladeX.ind \
	--consensus \
	--snp $f_snp \
	--groups \
		Wisent,Steppe,Aurochs,Bubalus_bubalis \
		CladeX,Steppe,Aurochs,Bubalus_bubalis \
	| tail -n +2 >> aabb-hypergeometric-cns.tsv

./aabb.py \
	--geno $f_geno \
	--ind Wisent_vs_CladeX.ind \
	--consensus \
	--snp $f_snp \
	--groups \
		Wisent,Aurochs,Steppe,Bubalus_bubalis \
		CladeX,Aurochs,Steppe,Bubalus_bubalis \
	| tail -n +2 >> aabb-hypergeometric-cns.tsv

###
exit





./aabb.py \
	--geno $f_geno \
	--ind $f_ind \
	--groups AncientWisent_4093A,Steppe_A3133,Aurochs_CPC98,Bubalus_bubalis \
		CladeX_6A,Steppe_A3133,Aurochs_CPC98,Bubalus_bubalis \
		HistoricalWisent_13947A,Steppe_A3133,Aurochs_CPC98,Bubalus_bubalis \
		ModernWisent,Steppe_A3133,Aurochs_CPC98,Bubalus_bubalis \
	> aabb-SteppeA3133-cns.tsv

./aabb.py \
	--geno $f_geno \
	--ind $f_ind \
	--groups AncientWisent_4093A,Steppe_A875,Aurochs_CPC98,Bubalus_bubalis \
		CladeX_6A,Steppe_A875,Aurochs_CPC98,Bubalus_bubalis \
		HistoricalWisent_13947A,Steppe_A875,Aurochs_CPC98,Bubalus_bubalis \
		ModernWisent,Steppe_A875,Aurochs_CPC98,Bubalus_bubalis \
	> aabb-SteppeA875-cns.tsv

./aabb.py \
	--geno $f_geno \
	--ind $f_ind \
	--groups AncientWisent,Aurochs,Steppe,Bubalus_bubalis \
		CladeX,Aurochs,Steppe,Bubalus_bubalis \
		HistoricalWisent,Aurochs,Steppe,Bubalus_bubalis \
		ModernWisent,Aurochs,Steppe,Bubalus_bubalis \
	> aabb-auroch.tsv

./aabb.py \
	--geno $f_geno \
	--ind $f_ind \
	--groups ModernWisent,AmericanBison,HistoricalWisent,Bubalus_bubalis \
		CladeX,AmericanBison,HistoricalWisent,Bubalus_bubalis \
		AncientWisent,AmericanBison,HistoricalWisent,Bubalus_bubalis \
	> aabb-american.tsv

