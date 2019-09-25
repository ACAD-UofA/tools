#!/bin/sh

module load grg-utils python/4.8.0/2.7.5 gnu/4.9.2

#./coverage -w121 /localscratch/jsoubrier/Bison_Genomes/CowRef/vcf/Bos_and_Bison.bcf.gz
#./coverage -w121 /localscratch/grg/bovid_genomes/bovids.bcf.gz

for hist in *.121.hist; do
	sample=${hist%.121.hist}
	plot_coverage.py \
		-w 121 \
		--title "$sample WindowSize=121" \
		$hist \
		$sample.cov121.pdf &
done

wait

