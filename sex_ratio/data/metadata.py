#!/usr/bin/env python

from __future__ import print_function
import sys
import csv
import re

def parse_sex(fn, skip=1):
    d = {}
    with open(fn) as f:
        while skip > 0:
            next(f)
            skip -= 1
        for line in f:
            fields = line.split()
            sample = fields[0]
            sex = fields[2]
            Nx, Na, Lx, La = map(int, fields[3:7])
            #if sex in "MF":
            d[sample] = (sex, Nx, Na, Lx, La)
    return d

def parse_elevation(fn):
    ll = {}
    with open(fn) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            ll[(row["lat"], row["lon"])] = row["elevation"]
    return ll

def parse_metadata_Bison_SHTGN(fn):

    with open(fn) as f:
        for row in csv.DictReader(f):#, delimiter="\t"):
            sample = row["SampleID"].split("_")[0]
            endog = row.get("Endogenous", row.get("ShotgunEndo"))
            clon = row.get("Clonal", row.get("ShotgunClonality"))

            yield sample, endog, clon

def parse_metadata_Bison(fn):

    # Site Type
    re_st_human = re.compile(r"(^|\s|[(-_.,])([Ss]tampede [Ss]ite|[Hh]unting)([)-_.,]|\s|$)")
    re_st_cave = re.compile(r"(^|\s|[(-_.,])([Cc]ave)([)-_.,]|\s|$)")

    # Sample Material
    re_sc_astr = re.compile(r"(^|\s|[(-_.,])([Aa]strag(alus)?|[Aa]stralagus|[Tt]alus)([)-_.,]|\s|$)")
    re_sc_foot = re.compile(r"(^|\s|[(-_.,])([Ff]oot|[Mm]etapod(ial)?|[Mm]etacar(pal|pus)?|[Mm]etatar(sal|sus)?|mx|mt|mc|meta(-)?x)([)-_.,]|\s|$)")
    re_sc_leg = re.compile(r"(^|\s|[(-_.,])([Hh]um(er(o)?us)?|[Rr]ad(ius)?|[Tt]ib(ia)?|[Ff]emur|[Uu]lna|[Ll]eg)([)-_.,]|\s|$)")
    re_sc_othercrania = re.compile(r"(^|\s|[(-_.,])([Bb]ulla|[Cc]rania|[Ss]kull|[Mm]andible|[Mm]axilla([re])?)([)-_.,]|\s|$)")
    re_sc_petrous = re.compile(r"(^|\s|[(-_.,])([Pp]etr(ous)?|[Pp]etrosal)([)-_.,]|\s|$)")
    re_sc_horn = re.compile(r"(^|\s|[(-_.,])([Hh]orn(core)?(s?))([)-_.,]|\s|$)")
    re_sc_tooth = re.compile(r"(^|\s|[(-_.,])([Tt]ooth)([)-_.,]|\s|$)")
    re_sc_vert = re.compile(r"(^|\s|[(-_.,])([Vv]ert(ebra(e)?)?|[Ss]pine|[Aa]tlas)([)-_.,]|\s|$)")
    re_sc_flat = re.compile(r"(^|\s|[(-_.,])([Pp]elvi[sc]|[Ss]capula)([)-_.,]|\s|$)")

    with open(fn) as f:
        for row in csv.DictReader(f):
            sample = row.get("ACAD Sample Number", row.get("ACAD Sample #")).strip()
            if sample[0] != "A":
                sample = "A" + sample

            loc = row.get("Locality")
            sloc = row.get("Specific Location", row.get("Englishized"))
            country = row.get("Country")
            scat1 = row.get("Sample Description", row.get("Sample type"))
            scat2 = row.get("Sample Details", row.get("Sample description"))

            museum = row.get("Museum Name", row.get("Museum"))
            museum_acc = row.get("Museum Accession Number", row.get("Museum #"))

            if re_st_cave.search(loc) or re_st_cave.search(sloc) or re_st_cave.search(scat2):
                site_type = "Cave"
            elif re_st_human.search(loc) or re_st_human.search(sloc) or re_st_human.search(scat2):
                site_type = "Anthropogenic"
                # skip these
                continue
            else:
                site_type = "Sediment"


            if re_sc_astr.search(scat1) or re_sc_astr.search(scat2):
                material = "Astragalus"
            elif re_sc_petrous.search(scat1) or re_sc_petrous.search(scat2):
                material = "Petrosal"
            elif re_sc_othercrania.search(scat1) or re_sc_othercrania.search(scat2):
                material = "OtherCrania"
            elif re_sc_leg.search(scat1) or re_sc_leg.search(scat2):
                material = "Leg"
            elif re_sc_foot.search(scat1) or re_sc_foot.search(scat2):
                material = "Foot"
            elif re_sc_horn.search(scat1) or re_sc_horn.search(scat2):
                material = "Horn"
            elif re_sc_tooth.search(scat1) or re_sc_tooth.search(scat2):
                material = "Tooth"
            elif re_sc_vert.search(scat1) or re_sc_vert.search(scat2):
                material = "Vertebra"
            elif re_sc_flat.search(scat1) or re_sc_flat.search(scat2):
                material = "FlatBone"
            else:
                material = None

            if material in {"Petrosal", "OtherCrania", "Horn", "Tooth"}:
                material2 = "Crania"
            elif material is not None:
                material2 = "Postcrania"
            else:
                material2 = None


            try:
                lat = float(row.get("Latitude", row.get("Latitude (N)")))
            except ValueError:
                lat = None
            try:
                lon = float(row.get("Longitude", row.get("Longitude (E)")))
            except ValueError:
                lon = None 

            #howloc = row.get("How location was obtained")

            age = row.get("OxCal4.2 Int13 Mean calBP", row.get("Carbon date (calibrated-mean)"))
            if age.endswith("calBP"):
                age = age[:-len("calBP")]
            elif age.endswith("BP"):
                age = age[:-len("BP")]
            if not age:
                age = None

            try:
                age = int(age)
            except (ValueError, TypeError):
                age = None

            endog = row.get("Shotgun Screening %Endogenous")
            try:
                endog = float(endog)
            except (ValueError, TypeError):
                endog = None

            clon = row.get("Shotgun Screening %Clonality")

            yield (sample, site_type, country, loc, sloc, material, material2, lat, lon, age, endog, museum, museum_acc)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("usage: {} sex.txt metadata.csv metadata2.csv latlonalt.txt > out.csv".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    sexdict = parse_sex(sys.argv[1])
    metadata = parse_metadata_Bison(sys.argv[2])
    shtgn = {}
    for sample, endog, clon in parse_metadata_Bison_SHTGN(sys.argv[3]):
        shtgn[sample] = (endog, clon)

    ll = parse_elevation(sys.argv[4])

    samples = set()

    print("sample", "sex", "site_type", "material", "material2", "lat", "lon", "alt", "age", "endog", "Nx", "Na", "Lx", "La",
            "Country", "Locality", "Specific Location", "Museum", "Museum Accession", sep="\t")
    for sample, site_type, country, loc, sloc, material, material2, lat, lon, age, endog, museum, museum_acc in metadata:
        sexfields = sexdict.get(sample, None)
        if sexfields is None:
            continue
        sex, Nx, Na, Lx, La = sexfields
        if sample in samples:
            continue
        samples.add(sample)
        if endog is None:
            endog, clon = shtgn.get(sample, (None, None))
        alt = ll.get((str(lat), str(lon)))

        if lon is not None and lon < -15:
            # Negative values are in Beringia and North America,
            # which are closer to Russia than to England, so make
            # the numbers reflect this.
            lon += 360

        print(sample, sex, site_type, material, material2, lat, lon, alt, age, endog, Nx, Na, Lx, La, country, loc, sloc, museum, museum_acc, sep="\t")
