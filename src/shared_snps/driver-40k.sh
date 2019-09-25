#!/bin/sh

idir=/localscratch/grg/40k/AmericanBison/02_all_40k_eig/
pfx=`cat /localscratch/grg/40k/AmericanBison/prefix`
f_geno=${idir}/${pfx}.geno
f_ind=${idir}/${pfx}.ind
f_snp=${idir}/${pfx}.snp

./aabb.py \
	--geno $f_geno \
	--ind $f_ind \
	--snp $f_snp \
	--groups \
		European,Steppe,Aurochs,Bubalus_bubalis \
		SteppeMN,Steppe,Aurochs,Bubalus_bubalis \
	> hypergeometric-40k.tsv

./aabb.py \
	--geno $f_geno \
	--ind $f_ind \
	--snp $f_snp \
	--groups \
		European,Aurochs,Steppe,Bubalus_bubalis \
		SteppeMN,Aurochs,Steppe,Bubalus_bubalis \
		Wood,Aurochs,Steppe,Bubalus_bubalis \
		Plains,Aurochs,Steppe,Bubalus_bubalis \
	>> hypergeometric-40k.tsv

./aabb.py \
	--geno $f_geno \
	--ind $f_ind \
	--snp $f_snp \
	--groups \
		European,Bubalus_bubalis,Steppe,Ovis_aries \
		SteppeMN,Bubalus_bubalis,Steppe,Ovis_aries \
		Wood,Bubalus_bubalis,Steppe,Ovis_aries \
		Plains,Bubalus_bubalis,Steppe,Ovis_aries \
	>> hypergeometric-40k.tsv
