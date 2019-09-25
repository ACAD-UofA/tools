#!/usr/bin/env python

from __future__ import print_function
import sys
import collections

def parse_snps_from_vcf(fn):
    if fn.endswith(".gz"):
        import gzip
        vopen = gzip.open
    else:
        vopen = open

    snpset = set()

    with vopen(fn) as f:
        for line in f:
            if line[0] == "#":
                continue
            fields = line.split("\t")
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            snpset.add((fields[0], fields[1], ref, alt))

    return snpset

samples = """
A16121_Bpriscus_Yukon_14792
A16171_Bpriscus_Yukon_40999
A1_BisonX_UralsRasik1_14874
A2494_Bprimigenius_RussiaYenseiRiver_14082
A3020_Bpriscus_Minnesota_7642
A3133_Bpriscus_Yukon_30612
A860_Bpriscus_Alaska_33154
A875_Bpriscus_Siberia_50k
BBO3569_Bbonasus_Bialowieza_0
BBO3574_Bbonasus_Bialowieza_0
BB_20087_Bbison_bison_0
Q229_Bbonasus_0
SW18_Bbison_athabascae_0
""".strip().split("\n")

samples = ["Ovis_aries_0"]

snpdd = collections.defaultdict(lambda:collections.defaultdict(set))
chrlist = map(str, range(1,26+1))

for s in samples:
    for c in chrlist:
        #vcf_fn = "stage5/" + s + ".snps-no-indels." + c + ".vcf.gz"
        vcf_fn = "02_depth/" + s + ".snpwin60." + c + ".vcf.gz"
        snpdd[s][c] = parse_snps_from_vcf(vcf_fn)

        snpdd["union"][c] |= snpdd[s][c]
        snpdd[s]["total"] |= snpdd[s][c]

    snpdd["union"]["total"] |= snpdd[s]["total"]

samples.append("union")
chrlist.append("total")

print("sample", *chrlist, sep="\t")
for s in samples:
    print(s, *[len(snpdd[s][c]) for c in chrlist], sep="\t")


snptypes = collections.Counter()
for _, _, ref, alt in snpdd["union"]["total"]:
    snptypes[(ref, alt)] += 1

tiset = set([("C","T"), ("T","C"), ("G","A"), ("A","G")])
ti = 0
tv = 0
print("\n\nref\talt\tcount")
for ref in "ACGT":
    for alt in "ACGT":
        if ref==alt:
            continue

        n = snptypes[(ref,alt)]
        if (ref,alt) in tiset:
            ti += n
        else:
            tv += n

        print(ref, alt, n, sep="\t")

print("\n\n")
print("Ti", ti, sep="\t")
print("Tv", tv, sep="\t")
