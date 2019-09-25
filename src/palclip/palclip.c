/*
 * Copyright (c) 2015 Graham Gower <graham.gower@gmail.com>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <getopt.h>
#include <stdint.h>
#include <zlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread);

void
usage(char *argv0)
{
	fprintf(stderr, "usage: %s [-m minpal] in.fast{a,q}\n", argv0);
	exit(1);
}

/*
 * Write out a single sequence.
 */
void
kseq_write(kseq_t *ks, FILE *fs)
{
	if (ks->qual.l == 0) {
		if (ks->comment.l)
			fprintf(fs, ">%s %s\n%s\n", ks->name.s, ks->comment.s, ks->seq.s);
		else
			fprintf(fs, ">%s\n%s\n", ks->name.s, ks->seq.s);
	} else {
		if (ks->comment.l)
			fprintf(fs, "@%s %s\n%s\n+\n%s\n", ks->name.s, ks->comment.s, ks->seq.s, ks->qual.s);
		else
			fprintf(fs, "@%s\n%s\n+\n%s\n", ks->name.s, ks->seq.s, ks->qual.s);
	}
}

static char rvmap[] = {['A'] = 'T', ['C'] = 'G', ['G'] = 'C', ['T'] = 'A'};

int
main(int argc, char **argv)
{
	char *f_in;
	int opt;
	unsigned int minpal = 4;

	gzFile fp;
	kseq_t *seq;
	int len;
	int i;
	uint64_t palcount = 0;

	while ((opt = getopt(argc, argv, "m:")) != -1) {
		switch (opt) {
			case 'm':
				errno = 0;
				minpal = strtoul(optarg, NULL, 0);
				if (errno || minpal > 1000)
					usage(argv[0]);
				break;
			default:
				usage(argv[0]);
		}
	}

	if (argc-optind != 1)
		usage(argv[0]);

	f_in = argv[optind];

	fp = gzopen(f_in, "r");
	seq = kseq_init(fp);
	while ((len = kseq_read(seq)) >= 0) {
		for (i=0; i<len/2; i++) {
			if (seq->seq.s[i] != rvmap[(int)seq->seq.s[len-i]])
				break;
		}

		if (i == len/2) {
			// whole sequence is a palindrome
			fprintf(stderr, "%s is a complete palindrome.\n", seq->name.s);
		}
		else if (i >= minpal) {
			seq->seq.l = len-i;
			if (seq->qual.l != 0)
				seq->qual.l = seq->seq.l;
			palcount++;
		}

		kseq_write(seq, stdout);
	}

	if (len != -1)
		fprintf(stderr, "%s: unexpected end of file.\n", f_in);

	fprintf(stderr, "found %jd palindromes of at least %dbp.\n", (intmax_t)palcount, minpal);

	kseq_destroy(seq);
	gzclose(fp);

	return 0;
}
