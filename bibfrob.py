#!/usr/bin/env python3

import sys
import re

# XXX: There is a biblib in pypi, but it is not the one I've used.
# Use this one instead:
# https://github.com/aclements/biblib.git
import biblib.bib

def parse_tex(fn):
    """
    Read LaTeX file to extract citation keys.

    This is a bit nasty, as it reads in the whole file in one go.
    """
    with open(fn) as f:
        latex = f.read()

    matches = re.findall(r"\\cite[pt](?:\[.*?\])?(?:\[.*?\])?\{(.*?)\}", latex, flags=re.DOTALL)
    latex = None

    if matches is None:
        print("{}: no citations".format(fn), file=sys.stderr)
        exit(1)

    cite_keys = set()
    for m in matches:
        for cite in re.split(r",|\n", m):
            for i,c in enumerate(cite):
                if c == '%':
                    cite = cite[:i]
                    break
            cite = cite.strip()
            if not cite:
                continue
            cite_keys.add(cite)

    return cite_keys

def parse_nlmcat(fn):
    """
    Parse NLM catalog.

    Downloaded from:
    https://www.ncbi.nlm.nih.gov/nlmcatalog?term=currentlyindexed
    Select "Send to:" -> "File" -> Format="Full (text)".
    """

    abbrvdict = {}
    abbrv = journal = None
    num_re = re.compile(r"^([0-9]+\. )")

    with open(fn) as f:
        for line in f:
            line = line.rstrip()

            if not line:
                if None in (journal, abbrv):
                    continue
                if journal in abbrvdict:
                    #raise Exception("abbrvdict[{}]={}, but want to add abbrvdict[{}]={}".format(x, abbrvdict[x], x, abbrv))
                    continue
                abbrvdict[journal.lower()] = abbrv
                j2 = journal.replace(".", "")
                if j2 != journal:
                    abbrvdict[j2.lower()] = abbrv
                abbrv = journal = None
                continue

            m = num_re.match(line)
            if m is not None:
                # skip entry number
                line = line[len(m.group(1)):]

            if line.startswith("Title Abbreviation: "):
                abbrv = line[len("Title Abbreviation: "):]
            elif line.startswith("Title(s): "):
                journal = line[len("Title(s): "):]

    # other variants
    abbrvdict["PNAS".lower()] = "Proc Natl Acad Sci U S A"
    abbrvdict["Proceedings of the National Academy of Sciences".lower()] = "Proc Natl Acad Sci U S A"
    abbrvdict["Philosophical Transactions of the Royal Society of London B: Biological Sciences".lower()] = "Philos Trans R Soc Lond B Biol Sci"
    abbrvdict["Journal of computational biology".lower()] = "J Comput Biol"

    return abbrvdict

def parse_bib(fn):
    try:
        with open(fn) as f:
            p = biblib.bib.Parser().parse(f, log_fp=sys.stderr)
            db = p.get_entries()
        db = biblib.bib.resolve_crossrefs(db)
    except biblib.messages.InputError:
        sys.exit(1)

    return db

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="frob bib file")
    parser.add_argument("--tex", metavar="citelist.tex",
                                    help="extract cite keys from this file")
    parser.add_argument("--nlm", metavar="nlmcatalog.txt",
                                    help="NLM catalog of journals")
    parser.add_argument("-v", "--verbose", action="store_true", default=False, help="increase verbosity")
    parser.add_argument("bibfn", metavar="file.bib", help="input file")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    db = parse_bib(args.bibfn)

    if args.tex:
        cite_keys = parse_tex(args.tex)
    else:
        cite_keys = db.keys()

    def author_year(s):
        fields = s.split("_")
        return (fields[0], fields[-1], fields[1:-1])
    cite_keys = sorted(cite_keys, key=author_year)

    if args.nlm:
        abbrvdict = parse_nlmcat(args.nlm)
    else:
        abbrvdict = {}

    unwanted_fields = ("abstract", "copyright", "keywords", "file", "pmid", "pmcid", "issn")
    wanted_fields = ("doi", "volume", "number", "pages")

    recoverer = biblib.messages.InputErrorRecoverer()
    with recoverer:
        for c in cite_keys:
            ent = db.get(c)
            if ent is None:
                print("Warn: citation {} not found".format(c), file=sys.stderr)
                continue
            for u in unwanted_fields:
                if u in ent:
                    del ent[u]

            if ent.typ == "article":
                journal = ent.get("journal")
                if journal is None:
                    print("Warn: {}: missing `journal' field.".format(ent.key), file=sys.stderr)
                else:
                    journal = journal.replace(".", "")
                    journal = journal.replace("\&", "&")
                    journal = abbrvdict.get(journal.lower(), journal)
                    ent["journal"] = journal
                    #print(journal)

                    if args.verbose:
                       for w in wanted_fields:
                            if journal == "bioRxiv" and w in ("volume", "number"):
                                continue
                            if journal == "Sci Rep" and w in ("number"):
                                continue
                            if ent.get(w) is None:
                                print("Warn: {}: missing `{}' field.".format(ent.key, w), file=sys.stderr)

            print(ent.to_bib())
            print()
    recoverer.reraise()
