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
#include <stdint.h>
#include <unistd.h>
#include "kseq.h"
#include "khash.h"

KSEQ_INIT(int, read);

#define max(a,b) ((a) > (b) ? (a) : (b))

typedef struct {
	char *name, *comment, *qual;
} read_aux_t;
KHASH_MAP_INIT_STR(str, read_aux_t);

/*
 * Write out a single sequence.
 */
void
seq_write(const char *seq, const read_aux_t *r, FILE *fs)
{
	if (r->qual) {
		// fastq
		if (r->comment)
			fprintf(fs, "@%s %s\n%s\n+\n%s\n", r->name, r->comment, seq, r->qual);
		else
			fprintf(fs, "@%s\n%s\n+\n%s\n", r->name, seq, r->qual);
	} else {
		// fasta
		if (r->comment)
			fprintf(fs, ">%s %s\n%s\n", r->name, r->comment, seq);
		else
			fprintf(fs, ">%s\n%s\n", r->name, seq);
	}
}


/*
 * In place reverse complement of s.
 */
void
revcomp(char *s, size_t len)
{
	static char cmap[] = {['A']='T', ['C']='G', ['G']='C', ['T']='A',
				['N']='N', ['n']='N',
				['a']='t', ['c']='g', ['g']='c', ['t']='a'};
	int i, j;
	int tmp;

	for (i=0, j=len-1; i<len/2; i++, j--) {
		tmp = cmap[(int)s[i]];
		s[i] = cmap[(int)s[j]];
		s[j] = tmp;
	}

	if (len % 2 == 1)
		s[i] = cmap[(int)s[i]];
}

void
usage(char *argv0)
{
	fprintf(stderr, "usage: %s < in.fast{a,q} > out.fast{a,q}\n", argv0);
	exit(1);
}

int
main(int argc, char **argv)
{
	kseq_t *seq;
	uint64_t dupes = 0, unique = 0, total = 0;
	khash_t(str) *h;
	int ret;
	int len;
	khint_t k;

	if (argc != 1)
		usage(argv[0]);

	seq = kseq_init(0);
	h = kh_init(str);

	while ((len = kseq_read(seq)) >= 0) {
		if (len == 0)
			continue;

		total++;

		// check forward strand
		if ((k = kh_get(str, h, seq->seq.s)) != kh_end(h)) {
			// found duplicate
			read_aux_t *r = &kh_val(h, k);
			if (seq->qual.l && r->qual) {
				int i;
				// take the highest quality bases
				for (i=0; i<seq->qual.l; i++)
					r->qual[i] = max(r->qual[i], seq->qual.s[i]);
			}
			dupes++;
			continue;
		}

		/*
		// check reverse strand
		revcomp(seq->seq.s, seq->seq.l);
		if ((k = kh_get(str, h, seq->seq.s)) != kh_end(h)) {
			// found duplicate
			read_aux_t *r = &kh_val(h, k);
			if (seq->qual.l && r->qual) {
				int i, j;
				// take the highest quality bases
				for (i=0, j=seq->qual.l-1; i<seq->qual.l; i++, j--)
					r->qual[i] = max(r->qual[i], seq->qual.s[j]);
			}
			dupes++;
			continue;
		}
		*/

		k = kh_put(str, h, seq->seq.s, &ret);
		if (ret > 0) {
			read_aux_t *r = &kh_val(h, k);
			r->name = strdup(seq->name.s);
			r->comment = seq->comment.l ? strdup(seq->comment.s) : NULL;
			r->qual = seq->qual.l ? strdup(seq->qual.s) : NULL;
			kh_key(h, k) = strdup(seq->seq.s);
		}
	}

	if (len != -1)
		fprintf(stderr, "Unexpected end of file.\n");

	// write out unique sequences
	for (k=kh_begin(h); k!=kh_end(h); k++) {
		if (kh_exist(h, k)) {
			char *seq = (char *)kh_key(h, k);
			read_aux_t *r = &kh_val(h, k);
			seq_write(seq, r, stdout);
			free(seq);
			free(r->name);
			if (r->comment)
				free(r->comment);
			if (r->qual)
				free(r->qual);
			unique++;
		}
	}

	kh_destroy(str, h);
	kseq_destroy(seq);

	fprintf(stderr, "total=%jd, unique=%jd, duplicates=%jd.\n",
			(intmax_t)total, (intmax_t)unique, (intmax_t)dupes);

	return 0;
}
