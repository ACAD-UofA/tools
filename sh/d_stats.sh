#!/bin/sh
# Calculate 4 population D(X, Y, Z, O) statistic for populations in $popfile.

tmpdir=/localscratch/grg/tmp
parfile=`mktemp $tmpdir/qpDstat.XXXXXX.par`
popfile=`mktemp $tmpdir/qpDstat.XXXXXX.popfilename`

module load admixtools

# contains 4 populations on each line <X> <Y> <Z> <O>
cat >$popfile << _EOF
European BisonX Steppe Ovis_aries
European Steppe BisonX Ovis_aries
BisonX Steppe European Ovis_aries

European Caucasicus Steppe Ovis_aries
European Steppe Caucasicus Ovis_aries
Caucasicus Steppe European Ovis_aries

BisonX Steppe Caucasicus Ovis_aries
Caucasicus Steppe BisonX Ovis_aries
BisonX Caucasicus Steppe Ovis_aries

Cow European Steppe Ovis_aries
Cow European SteppeMN Ovis_aries
Cow European Wood Ovis_aries
Cow European Wood_i Ovis_aries
Cow European Plains Ovis_aries
European Steppe Cow Ovis_aries
European SteppeMN Cow Ovis_aries
European Wood Cow Ovis_aries
European Wood_i Cow Ovis_aries
European Plains Cow Ovis_aries
Cow Steppe European Ovis_aries
Cow SteppeMN European Ovis_aries
Cow Wood European Ovis_aries
Cow Wood_i European Ovis_aries
Cow Plains European Ovis_aries

Cow Yak European Ovis_aries
Cow Yak Steppe Ovis_aries
Cow Yak SteppeMN Ovis_aries
Cow Yak Wood Ovis_aries
Cow Yak Wood_i Ovis_aries
Cow Yak Plains Ovis_aries
Cow Yak BisonX Ovis_aries
Cow Yak Caucasicus Ovis_aries
_EOF

pfx=9908_749Bovid_1ChimeraX_5SteppeMN_1Caucasicus_5Steppe_9BisonX
f_geno=${pfx}.geno
f_snp=${pfx}.snp
f_ind=${pfx}_grouplq.ind

cat >$parfile << _EOF
genotypename: $f_geno
snpname: $f_snp
indivname: $f_ind
popfilename: $popfile
_EOF

result2tsv() {
	echo -e "P1\tP2\tP3\tP4\tD\tZ\tBABA\tABBA\tN_SNPs"
	awk '$1=="result:" {
		for (i=2; i<=NF; i++) {
			if (i!=2)
				printf "\t";
			printf $i;
		}
		printf "\n";
	}'
}

qpDstat -p $parfile | result2tsv

rm $parfile $popfile
