#!/usr/bin/env python

import operator
import collections
import os.path

from scipy import stats

import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def parse_angsd_dstats(filename):
    data = []
    with open(filename) as f:
        next(f) # skip header
        for line in f:
            line = line.rstrip("\r\n")
            fields = line.split("\t")

            p1 = fields[0]
            p2 = fields[1]
            p3 = fields[2]
            #abba = int(fields[3])
            #baba = int(fields[4])
            d = float(fields[5])
            #jackEst = float(fields[6])
            #se = float(fields[7])
            z = float(fields[8])

            data.append((p1, p2, p3, d, z))
    return data

def parse_admixtools_dstats(filename):
    data = []
    with open(filename) as f:
        next(f) # skip header
        for line in f:
            line = line.rstrip("\r\n")
            fields = line.split("\t")

            p1 = fields[0]
            p2 = fields[1]
            p3 = fields[2]
            #p4 = fields[3]
            d = float(fields[4])
            z = float(fields[5])
            #baba = int(fields[6])
            #abba = int(fields[7])
            #n_sites = int(fields[8])

            # AdmixTools uses BABA-ABBA convention, which is opposite to that
            # of ANGSD, so flip the sign here so they match.
            d = -d
            z = -z

            data.append((p1, p2, p3, d, z))
    return data

def parse_dstats(filename):
    with open(filename) as f:
        header = next(f).rstrip("\r\n").split("\t")
    if header[0] == "H1":
        return parse_angsd_dstats(filename)
    elif header[0] == "P1":
        return parse_admixtools_dstats(filename)
    else:
        raise Exception("unknown input format for {}".format(filename))

def clean_dstats(data, dset):
    """
    Remove p1,p2,p3 combos not in dset; ensure correct p1,p2 orientation; sort.
    """

    data2 = {}
    for p1, p2, p3, d, z in data:
        if (p1,p2,p3) in dset:
            data2[(p1,p2,p3)] = (d,z)
        elif (p2,p1,p3) in dset:
            data2[(p2,p1,p3)] = (-d,-z)

    klist = sorted(data2.keys(), key=operator.itemgetter(0))
    klist.sort(key=operator.itemgetter(1))
    klist.sort(key=operator.itemgetter(2))
    return [(k[0],k[1],k[2],data2[k][0], data2[k][1]) for k in klist]

def get_unique_set(data, priordset=set()):
    """
    Remove identical stats, as D(p1,p2,p3,O) == -D(p2,p1,p3,O)

    If given a choice, use the orientation from priordset.
    """
    dset = set()
    for p1, p2, p3, _, _ in data:
        if not ((p1,p2,p3) in dset or (p2,p1,p3) in dset):
            if (p2,p1,p3) in priordset:
                dset.add((p2,p1,p3))
            else:
                dset.add((p1,p2,p3))
    return dset


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Compare two D-stats runs")
    parser.add_argument("tsv1", help="D-stats tsv file 1")
    parser.add_argument("tsv2", help="D-stats tsv file 2")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    d1 = parse_dstats(args.tsv1)
    d2 = parse_dstats(args.tsv2)

    dset1 = get_unique_set(d1)
    dset2 = get_unique_set(d2, dset1)

    dset = dset1 & dset2
    print(len(dset))

    d1 = clean_dstats(d1, dset)
    d2 = clean_dstats(d2, dset)

    outliers = collections.Counter()

    out_x1 = []
    out_x2 = []
    in_x1 = []
    in_x2 = []

    with open("/tmp/outliers.txt", "w") as f_outlier:
        for i, (dd1, dd2) in enumerate(zip(d1, d2)):
            if dd1[0:3] != dd2[0:3]:
                raise Exception("mismatch for entry {}:\n[1]=`{}`\n[2]=`{}`".format(i, dd1, dd2))

            if (((dd1[3] < 0 and dd2[3] > 0) or (dd1[3] > 0 and dd2[3] < 0))
                and (abs(dd1[4]) > 3 and abs(dd2[4]) > 3)):
                    f_outlier.write("[1]=`{}`\n[2]=`{}`\n".format(dd1, dd2))
                    for p in dd1[0:3]:
                        outliers[p] += 1
                    out_x1.append(dd1[3])
                    out_x2.append(dd2[3])
            else:
                in_x1.append(dd1[3])
                in_x2.append(dd2[3])

    print("{} outliers\n{}".format(sum(outliers.values())/3, outliers.most_common(10)))

    x1 = [dd[3] for dd in d1]
    x2 = [dd[3] for dd in d2]

    rho, pval = stats.spearmanr(x1, x2)
    print("spearman: rho={}, p={}".format(rho, pval))


    pdf = PdfPages("/tmp/test.pdf")
    #fig_w, fig_h = plt.figaspect(9.0/16.0)
    fig_w, fig_h = plt.figaspect(3.0/4.0)
    fig1 = plt.figure(figsize=(fig_w, fig_h))
    ax1 = fig1.add_subplot(111)

    ax1.scatter(in_x1, in_x2, color="red")
    ax1.scatter(out_x1, out_x2, color="red", edgecolor="blue")
    ax1.plot((-1,1), (-1,1), color="black", lw=2)
    ax1.set_xlim(-1,1)
    ax1.set_ylim(-1,1)
    ax1.set_xlabel(os.path.basename(args.tsv1))
    ax1.set_ylabel(os.path.basename(args.tsv2))

    ax1.annotate("spearman $\\rho={:.3f}$".format(rho),
            xy=(.99,-.99), ha='right', va='bottom')

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
