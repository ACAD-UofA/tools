#!/usr/bin/env python

from __future__ import print_function
import sys
import collections

def parse_mres(fn):
    m = collections.defaultdict(lambda:collections.defaultdict(list))
    with open(fn) as f:
        #data	model	Df	Deviance	Resid. Df	Resid. Dev	Pr(>Chi)
        next(f) # skip header
        for line in f:
            fields = line.rstrip().split("\t")
            data = fields[0]
            model = fields[1]
            df = fields[2]
            dev = "{:.3g}".format(float(fields[3]))
            resdf = fields[4]
            resdev = "{:.3g}".format(float(fields[5]))

            p = float(fields[6])
            if p < 0.05:
                #p = "\\multicolumn{{1}}{{Z{{.}}{{.}}{{-1}} }}{{{:.3G}}}".format(p)
                p = "\\bfseries {:.3g}".format(p)
            else:
                p = "{:.3g}".format(p)

            m[model][data].append((df, dev, resdf, resdev, p))
    return m

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: {} mres.txt".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    rows = ("site_type", "material", "material2", "age",
                "lat", "lon", "alt", "alpine", "endog", "GC", "readlen", "C2T")
    cols = ("bison", "bison (postcrania)", "bison (sediment)",
            "brown bears", "brown bears (Alps)", "brown bears (non Alps)",
            "mammoths")

    crename = {"bison (sediment)": "bison (non-cave)"}

    rrename = { "intercept-only": "Intercept-only",
                "site_type": "Cave/non-cave",
                "material": "Material1",
                "material2": "Material2",
                "age": "$^{{14}}$C age",
                "lat": "Latitude",
                "lon": "Longitude",
                "alt": "Altitude",
                "alpine": "Alps/non-Alps",
                "endog": "Endogenous DNA",
                "GC": "GC ratio",
                "readlen": "DNA fragment length",
                "C2T": "5' C$\\rightarrow$T"}

    rcat = set(["site_type", "material", "material2", "alpine"])

    m = parse_mres(sys.argv[1])

    # data	model	coeff	Estimate	Std. Error	z value	Pr(>|z|)
    #print("Est.", "S.E.", "$Z$", "$p$", sep=" & ", end="\\\\ ")
    #print("\\midrule\\\\")
    #data	model	Df	Deviance	Resid. Df	Resid. Dev	Pr(>Chi)

    for i, col in enumerate(cols,1):
        c = crename.get(col, col)
        for j, row in enumerate(rows,1):
            mrc = m[row][col]
            if len(mrc) == 0:
                continue

            if row == "intercept-only":
                sym = ""
            elif row in rcat:
                sym = "$^\dag$"
            else:
                sym = "$^\ddag$"

            r = sym + rrename[row]
            for df, dev, resdf, resdev, p in mrc:
                #if float(p) < 0.05:
                #    print("\\rowstyle{\\bfseries\\boldmath} ", end="")
                print(c, r, df, dev, resdf, resdev, p, sep=" & ", end="\\\\\n")

            if j != len(rows):
                print("\\noalign{\\smallskip}")

        if i != len(cols):
            print("\\noalign{\\bigskip}")

