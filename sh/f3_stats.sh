#!/bin/sh
# Calculate 3 population test f3(C; A, B) for populations in $popfile.

tmpdir=/localscratch/grg/tmp
parfile=`mktemp $tmpdir/qpDstat.XXXXXX.par`
popfile=`mktemp $tmpdir/qpDstat.XXXXXX.popfilename`

module load admixtools

# contains 3 populations on each line <Source1 (A)> <Source2 (B)> < Target (C)>
cat >$popfile << _EOF
Wood Cow Wood_i
Wood Wood_i Cow
Wood_i Cow Wood
Wood Plains Cow
Bos_Gaurus Cow Zebu
Yak European Wood
_EOF

pfx=39833_749Bovid
f_geno=${pfx}.geno
f_snp=${pfx}.snp
f_ind=${pfx}.ind

cat >$parfile << _EOF
genotypename: $f_geno
snpname: $f_snp
indivname: $f_ind
popfilename: $popfile
_EOF


result2tsv() {
	# The results have the following format - 
	# result:   Source1  Source2   Target f_3  std.err Z SNPs
	echo -e "Source1\tSource2\tTarget\tF3\tStd.Err\tZ\tN_SNPs"
	awk '$1=="result:" {
		for (i=2; i<=NF; i++) {
			if (i!=2)
				printf "\t";
			printf $i;
		}
		printf "\n";
	}'
}


qp3Pop -p $parfile | result2tsv

rm $parfile $popfile
