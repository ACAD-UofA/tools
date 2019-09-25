import sys
import swalign

scoring = swalign.NucleotideScoringMatrix()
sw = swalign.LocalAlignment(scoring)

def selfanneal(s):
    aln = sw.align(s, swalign.revcomp(s))
    return aln.score

if __name__ == "__main__":
    for s in sys.argv[1:]:
        print(selfanneal)
