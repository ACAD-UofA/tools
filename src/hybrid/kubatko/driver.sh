
steppes_all="
Steppe
Steppe_bx
Steppe_A16121
Steppe_A860
Steppe_A875
Steppe_A885
Steppe_A3133
"
steppes="
Steppe_A875+A3133
Steppe_A875
Steppe_A3133
"

wisent="
CladeX
6A
HistoricalWisent
AncientWisent
Wisent+CladeX
"

for s in $steppes; do
	for w in $wisent; do
		echo "### $s $w Aurochs"
		./meng 100 20000 bovids.vcf bovids.txt $s $w Aurochs
	done
done
