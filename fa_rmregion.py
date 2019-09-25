#!/usr/bin/env python

def parse_fa(filename):
    contig_header = None
    seq = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line[0] == ">":
                if contig_header:
                    yield contig_header, "".join(seq)
                contig_header = line
                seq = []
            else:
                seq.append(line)

    # last sequence
    yield contig_header, "".join(seq)

def print_seq(contig_header, seq, linelen=60):
    print(contig_header)
    i = 0
    while i<len(seq):
        print(seq[i:i+linelen])
        i += linelen

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="delete region [a, b) from a fasta file")
    parser.add_argument("fasta_file", help="fasta file")
    parser.add_argument("contig", help="name of contig or chromosome")
    parser.add_argument("a", type=int, help="start position (0-based coordinates)")
    parser.add_argument("b", type=int, help="end position (0-based coordinates)")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    if args.a >= args.b:
        raise Exception("a({}) >= b({}) doesn't make sense".format(args.a, args.b))

    for contig_header, seq in parse_fa(args.fasta_file):
        fields = contig_header.split(" ")
        contig = fields[0][1:]
        if contig == args.contig:
            seq = seq[:args.a] + seq[args.b:]
        print_seq(contig_header, seq)
