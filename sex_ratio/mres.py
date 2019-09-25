#!/usr/bin/env python

from __future__ import print_function
import sys
import collections

def parse_mres(fn):
    m = collections.defaultdict(lambda:collections.defaultdict(list))
    with open(fn) as f:
        # data	model	coeff	Estimate	Std. Error	z value	Pr(>|z|)
        next(f) # skip header
        for line in f:
            fields = line.rstrip().split("\t")
            data = fields[0]
            model = fields[1]
            coeff = fields[2]

            est, se, z = ("{:.3g}".format(float(x)) for x in fields[3:6])

            p = float(fields[6])
            if p < 0.05:
                #p = "\\multicolumn{{1}}{{Z{{.}}{{.}}{{-1}} }}{{{:.3G}}}".format(p)
                p = "\\bfseries {:.3g}".format(p)
            else:
                p = "{:.3g}".format(p)

            #yield model, data, coeff, est, se, z, p
            m[model][data].append((coeff, est, se, z, p))
    return m

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: {} mres.txt".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    rows = ("intercept-only", "site_type", "material", "material2", "age",
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
            for i, (coeff, est, se, z, p) in enumerate(mrc):
                #if float(p) < 0.05:
                #    print("\\rowstyle{\\bfseries\\boldmath} ", end="")
                if i == 1 and row not in rcat:
                    coeff = "slope"
                print(c, r, coeff, est, se, z, p, sep=" & ", end="\\\\\n")

            if j != len(rows):
                print("\\noalign{\\smallskip}")

        if i != len(cols):
            print("\\noalign{\\bigskip}")

