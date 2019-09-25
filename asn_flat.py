#!/usr/bin/env python
#
# Parse dbSNP `ASN.1 flat' files from here:
# ftp://ftp.ncbi.nih.gov/snp/organisms/cow_9913/ASN1_flat/

from __future__ import print_function
import sys
from gzopen import gzopen

def parse_asn(fn, assembly):

    def subsplit(lst, tok="="):
        sdict = dict()
        for elm in lst:
            elm = elm.strip()
            subfields = elm.split(tok)
            if len(subfields) == 1:
                sdict[elm] = True
            elif len(subfields) == 2:
                sdict[subfields[0]] = subfields[1]
        return sdict

    rs = ss = chrom = pos = ref = alt = None
    indel = False

    with gzopen(fn) as f:
        # skip header
        next(f); next(f); next(f)
        for line in f:
            line = line.strip()

            if line == "":
                # next entry
                if None not in (rs, ss, chrom, pos):
                    yield rs, ss, chrom, pos, ref, alt
                rs = ss = chrom = pos = ref = alt = None
                indel = False
                continue

            if indel:
                continue

            fields = line.split(" | ")

            if fields[0].startswith("rs"):
                subfields = subsplit(fields[1:])
                if "snp" in subfields:
                    rs = fields[0]
                else:
                    indel = True

            elif fields[0].startswith("ss"):
                if ss is None:
                    ss = [fields[0]]
                else:
                    ss.append(fields[0])

            elif fields[0] == "SNP":
                subfields = subsplit(fields[1:])
                alleles = subfields.get("alleles")
                if alleles is not None:
                    if len(alleles) == 5:
                        # biallelic string, e.g. 'A/G'
                        ref, alt = alleles[1], alleles[3]

            elif fields[0] == "CTG" and chrom is None:
                subfields = subsplit(fields[1:])
                if subfields.get("assembly") == assembly:
                    chrom = subfields.get("chr")
                    pos = int(subfields.get("chr-pos", -1))

    # last entry
    if None not in (rs, ss, chrom, pos):
        yield rs, ss, chrom, pos, ref, alt

def asn2dict(asn_flat_fn, assembly="Bos_taurus_UMD_3.1.1"):
    asn_dict = dict()
    for entry in parse_asn(asn_flat_fn, assembly):
        rs, ss, chrom, pos, ref, alt = entry
        asn_dict[rs] = entry
        for ss_id in ss:
            asn_dict[ss_id] = entry
    return asn_dict

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: {} ds_flat.gz".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    asn_flat_fn = sys.argv[1]

    for entry in parse_asn(asn_flat_fn, "Bos_taurus_UMD_3.1.1"):
        rs, ss, chrom, pos, ref, alt = entry
        print(rs, ss[0], chrom, pos, ref, alt, sep="\t")
        #for ss_id in ss:
        #    print(ss_id, rs, chrom, pos, ref, alt, sep="\t")
