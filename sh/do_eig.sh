#!/bin/sh

twstats=/localscratch/grg/twtable
prefix="9908_749Bovid_1ChimeraX_5SteppeMN_1Caucasicus_5Steppe_9BisonX"
popsets="
BisonX:BisonX_lq:ChimeraX:Caucasicus:Steppe:Steppe_16121:Steppe_lq:SteppeMN:SteppeMN_lq
BisonX:BisonX_lq:ChimeraX:Caucasicus:Steppe:Steppe_16121:Steppe_lq:SteppeMN:SteppeMN_lq:Ovis_aries
BisonX:BisonX_lq:ChimeraX:Caucasicus:European:Wood:Wood_i:Steppe:Steppe_16121:Steppe_lq:SteppeMN:SteppeMN_lq
BisonX:BisonX_lq:ChimeraX:Caucasicus:European:Wood:Wood_i:Steppe:Steppe_16121:Steppe_lq:SteppeMN:SteppeMN_lq:Yak
BisonX:BisonX_lq:ChimeraX:Caucasicus:European:Wood:Wood_i:Steppe:Steppe_16121:Steppe_lq:SteppeMN:SteppeMN_lq:Ovis_aries
BisonX:BisonX_lq:ChimeraX:Caucasicus:European:Wood:Wood_i:Steppe:Steppe_16121:Steppe_lq:SteppeMN:SteppeMN_lq:Ovis_aries:Yak
BisonX:BisonX_lq:ChimeraX:Caucasicus:European:Wood:Steppe:Steppe_16121:Steppe_lq:SteppeMN:SteppeMN_lq
BisonX:BisonX_lq:ChimeraX:Caucasicus:European:Wood:Steppe:Steppe_16121:Steppe_lq:SteppeMN:SteppeMN_lq:Yak
BisonX:BisonX_lq:ChimeraX:Caucasicus:European:Wood:Steppe:Steppe_16121:Steppe_lq:SteppeMN:SteppeMN_lq:Ovis_aries
BisonX:BisonX_lq:ChimeraX:Caucasicus:European:Wood:Steppe:Steppe_16121:Steppe_lq:SteppeMN:SteppeMN_lq:Ovis_aries:Yak
BisonX:Caucasicus:Steppe:SteppeMN
BisonX:Caucasicus:Steppe:SteppeMN:Ovis_aries
BisonX:Caucasicus:European:Wood:Steppe:SteppeMN
BisonX:Caucasicus:European:Wood:Steppe:SteppeMN:Yak
BisonX:Caucasicus:European:Wood:Steppe:SteppeMN:Ovis_aries
"

# eigensoft needs libgsl, gnuplot and ps2pdf from ghostscript
module load ghostscript gnuplot gsl eigensoft

#cd ./eig_6A_ChimeraX_with_Steppe

# ploteig uses a bad shebang path for perl
ploteig="perl -w `which ploteig`"

for popset in $popsets; do

	# create acronym from the popset
	popcode=`echo $popset | awk -F: '{for (i=1; i<=NF; i++) {printf "%s",substr($i, 0, 1)}}'`

	# put list of population groups into a file
	poplist="poplist_${popcode}.txt"
	echo $popset | awk -F: '{for (i=1; i<=NF; i++) print $i}' > $poplist

	oprefix="${prefix}_${popcode}"

	# calculate eigenvectors/eigenvalues for the first 4 principle components
	smartpca.perl \
		-i ${prefix}.geno \
		-a ${prefix}.snp \
		-b ${prefix}.ind \
		-o ${oprefix}.pca \
		-p ${oprefix}.plot \
		-e ${oprefix}.eval \
		-l ${oprefix}.log \
		-w $poplist \
		-y $poplist \
		-k 6 \
		-m 0
	# project onto PWE principle components with:
	# -w poplist_PWE.txt

	# Test for significance of principle components,
	# using Tracy-Widom theory.  See Patterson et al. 2006.
	twstats -t $twstats -i ${oprefix}.eval -o ${oprefix}.twstats

	# plot some principle component combinations
	for ev in "1:2" "3:4" "5:6"; do

		# Remove colon, for use in output filename.
		# Windows disallows colons in filenames and mangles the name.
		# Bash completion also has problems with colons.
		ev_no_colon=`echo $ev | sed 's/://'`
		
		$ploteig \
			-c $ev \
			-i ${oprefix}.pca.evec \
			-p $poplist \
			-o ${oprefix}_${ev_no_colon}.xtxt \
			-f \
			-x
	done
done
