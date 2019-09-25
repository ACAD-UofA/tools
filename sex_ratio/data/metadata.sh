#!/bin/sh

./metadata.py \
	AmBison.txt \
	LEF_Bison/AmericanBisonatACAD_10_01_17-1.csv \
	LEF_Bison/Bisonallalvl_overview_Box2017_1.SHTGN_stats.csv \
	latlonalt.txt \
	> am-bison.txt

./metadata.py \
	EuBison.txt \
	LEF_Bison/Bisonallalvl_overview_Box2017_1.csv \
	LEF_Bison/Bisonallalvl_overview_Box2017_1.SHTGN_stats.csv \
	latlonalt.txt \
	> eu-bison.txt

(cat am-bison.txt; tail -n +2 eu-bison.txt) > bison.txt
Rscript bison-data-frob.R
rm am-bison.txt eu-bison.txt bison.txt


awk -f metadata_bears.awk \
	< Alex_GLM_Bear/SexRatioBrownBear2.0_with_alt.csv \
	> meta.bears.txt
Rscript brownbear-data-frob.R
rm meta.bears.txt


Rscript mammoth-data-frob.R
