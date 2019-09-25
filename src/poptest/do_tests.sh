#!/bin/sh

vcf1=/localscratch/grg/H_sapiens/1Kgenomes/vcf/1KG_and_ASA_10_or_more_ASA_calls.vcf
vcf2=/localscratch/grg/H_sapiens/1Kgenomes/vcf/1KG_and_ASA_10_or_more_ASA_Europeans_calls.vcf
vcf3=/localscratch/grg/H_sapiens/1Kgenomes/vcf/1KG.chr_prefix.vcf
odir=out

HTSDIR=../htslib
BCFDIR=../bcftools
export LD_LIBRARY_PATH=$HTSDIR:$LD_LIBRARY_PATH
export BCFTOOLS_PLUGINS=./
#valgrind --leak-check=full --show-leak-kinds=all

body() {
    IFS= read -r header
    printf '%s\n' "$header"
    "$@"
}

mkdir -p $odir

tmpfile=`mktemp`

# Ancient Americans vs. PEL
$BCFDIR/bcftools +poptest \
	$vcf1 \
	-- \
	-p /localscratch/grg/H_sapiens/SNP_Capture_Europe_Americas_ACAD.txt.Americas \
	-p /localscratch/grg/H_sapiens/1Kgenomes/1Kgenomes.PEL.ind \
	| body sort -n -k 6 \
	> $tmpfile

	./fst2zfst_append.py -c 9 $tmpfile \
		| ./annotation_append.py /dev/stdin > $odir/Americans_vs_PEL.2p.tsv

# Ancient Americans vs. CEU
$BCFDIR/bcftools +poptest \
	$vcf1 \
	-- \
	-p /localscratch/grg/H_sapiens/SNP_Capture_Europe_Americas_ACAD.txt.Americas \
	-p /localscratch/grg/H_sapiens/1Kgenomes/1Kgenomes.CEU.ind \
	| body sort -n -k 6 \
	> $tmpfile

	./fst2zfst_append.py -c 9 $tmpfile \
		| ./annotation_append.py /dev/stdin > $odir/Americans_vs_CEU.2p.tsv

# Ancient Europeans vs. CEU
$BCFDIR/bcftools +poptest \
	$vcf2 \
	-- \
	-p /localscratch/grg/H_sapiens/SNP_Capture_Europe_Americas_ACAD.txt.Europe \
	-p /localscratch/grg/H_sapiens/1Kgenomes/1Kgenomes.CEU.ind \
	| body sort -n -k 6 \
	> $tmpfile

	./fst2zfst_append.py -c 9 $tmpfile | \
		./annotation_append.py /dev/stdin > $odir/Europeans_vs_CEU.2p.tsv

# Ancient Europeans vs. PEL
$BCFDIR/bcftools +poptest \
	$vcf2 \
	-- \
	-p /localscratch/grg/H_sapiens/SNP_Capture_Europe_Americas_ACAD.txt.Europe \
	-p /localscratch/grg/H_sapiens/1Kgenomes/1Kgenomes.PEL.ind \
	| body sort -n -k 6 \
	> $tmpfile

	./fst2zfst_append.py -c 9 $tmpfile | \
		./annotation_append.py /dev/stdin > $odir/Europeans_vs_PEL.2p.tsv

# CEU vs. PEL
$BCFDIR/bcftools +poptest \
	$vcf3 \
	-- \
	-p /localscratch/grg/H_sapiens/1Kgenomes/1Kgenomes.CEU.ind \
	-p /localscratch/grg/H_sapiens/1Kgenomes/1Kgenomes.PEL.ind \
	| body sort -n -k 6 \
	> $tmpfile

	./fst2zfst_append.py -c 9 $tmpfile | \
		./annotation_append.py /dev/stdin > $odir/CEU_vs_PEL.2p.tsv

# CEU vs. CLM
$BCFDIR/bcftools +poptest \
	$vcf3 \
	-- \
	-p /localscratch/grg/H_sapiens/1Kgenomes/1Kgenomes.CEU.ind \
	-p /localscratch/grg/H_sapiens/1Kgenomes/1Kgenomes.CLM.ind \
	| body sort -n -k 6 \
	> $tmpfile

	./fst2zfst_append.py -c 9 $tmpfile | \
		./annotation_append.py /dev/stdin > $odir/CEU_vs_CLM.2p.tsv

# CEU vs. MXL
$BCFDIR/bcftools +poptest \
	$vcf3 \
	-- \
	-p /localscratch/grg/H_sapiens/1Kgenomes/1Kgenomes.CEU.ind \
	-p /localscratch/grg/H_sapiens/1Kgenomes/1Kgenomes.MXL.ind \
	| body sort -n -k 6 \
	> $tmpfile

	./fst2zfst_append.py -c 9 $tmpfile | \
		./annotation_append.py /dev/stdin > $odir/CEU_vs_MXL.2p.tsv

rm -f $tmpfile
