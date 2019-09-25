#!/usr/bin/env python

from __future__ import print_function
import sys
import csv
from operator import itemgetter
from gzopen import gzopen

def parse_trusted(fn):
    rev = {"ADBR": True, "ARBD": False}
    with gzopen(fn) as f:
        for row in csv.reader(f):
            snp_num = int(row[0])
            chrom = row[1]
            if chrom == "30":
                chrom = "X"
            pos = int(row[2])
            ref = row[3]
            alt = row[4]
            key = row[5]
            yield snp_num, chrom, pos, ref, alt, rev[key]

def parse_dbSNP_dump(fn):
    with gzopen(fn) as f:
        for line in f:
            fields = line.split()
            rs = fields[0]
            ss = fields[1]
            chrom = fields[2]
            pos = int(fields[3])
            ref = fields[4]
            alt = fields[5]
            yield rs, ss, chrom, pos, ref, alt

def parse_snp_chip_txt(fn):
    data = {}
    lineno = 0
    with gzopen(fn) as f:
        lineno += 1
        hline = next(f).rstrip("\r\n")
        headers = hline.split("\t")[1:]
        for line in f:
            line = line.rstrip("\r\n")
            fields = line.split("\t")
            snp_num = int(fields[0])
            data[snp_num] = fields[1:]

            if len(data[snp_num]) != len(headers):
                raise ParseError("Error: {}: line {}: mismatch between number of headers ({}) and data array length ({}).".format(filename, lineno, len(headers), len(data[snp_num])))

    return headers, data

def parse_snp50_csv(fn):
    with gzopen(fn) as f:
        for row in csv.DictReader(f):
            snp_num = int(row["snp_number"])
            ss = row["ss_id"]
            rs = row["rs_id"]
            if rs:
                rs = "rs"+rs
            chrom = row["umd30_bta"]
            if chrom == "30":
                chrom = "X"
            pos = row["umd30_pos"]
            try:
                pos = int(pos)
            except ValueError:
                pass
            yield snp_num, ss, rs, chrom, pos

def sort_loci_by_coord(loci_list):
    def chr_cmp(a, b):
        if a == b:
            return 0
        try:
            if int(a) < int(b):
                return -1
            else:
                return 1
        except ValueError:
            return cmp(a, b)
    loci_list.sort(key=itemgetter(1))
    loci_list.sort(cmp=chr_cmp, key=itemgetter(0))

def eig_chr_map(chrom, nchrom=29):
    if chrom == "X":
        return nchrom+1
    elif chrom == "Y":
        return nchrom+2
    elif chrom in ("M","MT","Mt"):
        return 90
    else:
        try:
            chrom_int = int(chrom)
        except ValueError:
            raise RuntimeError("Unknown chromosome {}".format(chrom_str))
        return chrom_int

def write_indfile(bnames, oprefix):
    with open("{}.ind".format(oprefix), "w") as f:
        for h in bnames:
            hfields = h.split("_")

            # Come up with (semi) sensible group/sample names.
            if len(hfields) == 1:
                group = h
            else:
                if h[0] != "X":
                    group = "_".join(hfields[:-1])
                else:
                    if hfields[1][0] in ("A", "B"):
                        group = "_".join(hfields[2:])
                    else:
                        group = h
            if group.startswith("rep_"):
                group = group[4:]
            elif group.startswith("replig_"):
                group = group[7:]
            elif group.startswith("pool_"):
                group = group[5:]
            elif group.startswith("lig_col_"):
                group = group[8:]
            elif group.startswith("lig_eth_"):
                group = group[8:]
            elif group.startswith("3ul_"):
                group = group[4:]
            elif group.startswith("6ul_"):
                group = group[4:]
            elif group.startswith("10ul_"):
                group = group[5:]

            f.write("{}\tU\t{}\n".format(h, group))

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Convert SNPs from SNP chip .txt to EIGENSOFT format (with 2=REF/REF, 1=REF/ALT, 0=ALT/ALT, 9=unknown).")
    parser.add_argument("-o", "--oprefix", required=True, help="Output prefix for eigenstrat *.{geno,ind,snp} files")
    parser.add_argument("transposed_fn", help="SNP chip data (40843_749Bovid_Transposed.txt.gz)")
    parser.add_argument("snp50_csv_fn", help="50k SNP info (SNP50_complete_info.csv.gz)")
    parser.add_argument("trusted_fn", help="Trusted SNPs, with encoding direction (trusted_chip2seq_nov8_14.txt.gz)")
    parser.add_argument("dbSNP_dump", help="my dbSNP dump (rs_ss0_coords.txt.gz)")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    bnames, genotypes = parse_snp_chip_txt(args.transposed_fn)
    write_indfile(bnames, args.oprefix)

    #snp_num, chrom, pos, ref, alt, rev[key]
    trusted = {row[0]:row for row in parse_trusted(args.trusted_fn)}

    dbSNP_by_coord = dict()
    dbSNP_by_rs = dict()
    dbSNP_by_ss = dict()
    for row in parse_dbSNP_dump(args.dbSNP_dump):
        rs, ss, chrom, pos, ref, alt = row
        dbSNP_by_coord[(chrom,pos)] = row
        dbSNP_by_rs[rs] = row
        dbSNP_by_ss[ss] = row

    gt_map_fwd = {"0":"2", "1":"1", "2":"0", "?":"9"}
    gt_map_rev = {"0":"0", "1":"1", "2":"2", "?":"9"}
    observed = set()

    validated = []

    for row in parse_snp50_csv(args.snp50_csv_fn):
        snp_num, ss, rs, chrom, pos = row

        gt_str = genotypes.get(snp_num)
        if gt_str is None:
            continue

        tr_row = trusted.get(snp_num)
        if tr_row is None:
            continue
        _, _, _, tr_ref, tr_alt, rev = tr_row

        db1 = dbSNP_by_coord.get((chrom,pos))
        db2 = dbSNP_by_rs.get(rs)
        db3 = dbSNP_by_ss.get(ss)
        if db1 is not None:
            db = db1
        elif db2 is not None:
            db = db2
        elif db3 is not None:
            db = db3
        else:
            print("Warning: {}: no valid dbSNP entry".format(snp_num), row, file=sys.stderr)
            continue
        db_rs, db_ss, db_chrom, db_pos, db_ref, db_alt = db

        skip = 0

        for x in set([db1,db2,db3]):
            if x is not None:
                if x in observed:
                    print("Warning: {}: duplicate entry".format(snp_num), x, file=sys.stderr)
                    skip = 1
                observed.add(x)

        if db1 is not None and db2 is not None:
            if db1 is not db2:
                print("Warning: {}: db1!=db2: ({}) != ({})".format(snp_num, db1, db2), file=sys.stderr)
                skip = 1
        if db2 is not None and db3 is not None:
            if db2 is not db3:
                print("Warning: {}: db2!=db3: ({}) != ({})".format(snp_num, db2, db3), file=sys.stderr)
                skip = 1
        if db1 is not None and db3 is not None:
            if db1 is not db3:
                print("Warning: {}: db1!=db3: ({}) != ({})".format(snp_num, db1, db3), file=sys.stderr)
                skip = 1

        if skip:
            continue

        if rev:
            gt_map = gt_map_rev
        else:
            gt_map = gt_map_fwd

        gts = "".join(gt_map[gt] for gt in gt_str)
        entry = (eig_chr_map(db_chrom), db_pos, db_rs, tr_ref, tr_alt, gts)
        validated.append(entry)

    validated.sort(key=itemgetter(0,1))

    with open("{}.snp".format(args.oprefix), "w") as f_snp, \
         open("{}.geno".format(args.oprefix), "w") as f_geno:
        for chrom, pos, rs, ref, alt, gts in validated:
            print(rs, chrom, 0.0, pos, ref, alt, sep="\t", file=f_snp)
            #print(rs, chrom, 0.0, pos, sep="\t", file=f_snp)
            print(gts, file=f_geno)

