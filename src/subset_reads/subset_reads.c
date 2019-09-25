/*
 * Create subsets of sequences from a fastq, with decreasing sequence lengths.
 *
 * Copyright (c) 2016 Graham Gower <graham.gower@gmail.com>
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
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdint.h>
#include <errno.h>
#include <time.h>

#include "kseq.h"
#include "kmath.h"

KSEQ_INIT(int, read);


int
subset(char *fn, int nrep, const int *sizes, int n)
{
	kseq_t *seq;
	int fd;
	int i, j;
	int len;
	int ret = 0;
	int s_off, s_len;
	uint64_t seqno;
	static krand_t *kr;

	fd = open(fn, O_RDONLY);
	if (fd == -1) {
		fprintf(stderr, "Error: %s: %s\n", fn, strerror(errno));
		ret = -1;
		goto err0;
	}

	seq = kseq_init(fd);

	kr = kr_srand(time(NULL)*getpid());

	seqno = 0;
	while ((len = kseq_read(seq)) >= 0) {
		seqno++;

		if (seq->qual.l == 0) {
			fprintf(stderr, "Error: %s: missing quality scores.\n", fn);
			ret = -5;
			goto err1;
		}

		// original sequence
		printf("@%s\n%s\n+\n%s\n", seq->name.s, seq->seq.s, seq->qual.s);

		for (i=0; i<nrep; i++) {
			s_off = 0;
			s_len = seq->seq.l;
			for (j=0; j<n; j++) {
				if (s_len < sizes[j]) {
					fprintf(stderr, "Error: %s: line %jd: sequence length (%d) < size[%d]=%d\n",
							fn, (uintmax_t)(seqno*4+1), s_len, j, sizes[j]);
					ret = -6;
					goto err1;
				}
				s_off += kr_rand(kr) % (s_len - sizes[j]);
				s_len = sizes[j];

				printf("@%s:%d:%d\n", seq->name.s, i+1, s_len);
				fwrite(seq->seq.s+s_off, 1, s_len, stdout);
				fwrite("\n+\n", 1, 3, stdout);
				fwrite(seq->qual.s+s_off, 1, s_len, stdout);
				fwrite("\n", 1, 1, stdout);
			}
		}
	}

	if (len != -1) {
		fprintf(stderr, "Warning : %s: unexpected end of file.\n", fn);
	}

err1:
	kseq_destroy(seq);
	close(fd);
err0:
	return ret;
}

int
cmp_n(const void *p1, const void *p2)
{
	int n1 = *(intptr_t*)p1;
	int n2 = *(intptr_t*)p2;
	return n2-n1;
}

int
main(int argc, char **argv)
{
	int *sizes;
	int nrep;
	int n;
	int i;

	if (argc < 4) {
		fprintf(stderr, "usage: %s in.fastq NREP SZ1 [SZ2 ... SZx]\n", argv[0]);
		return -1;
	}

	errno = 0;
	nrep = strtoul(argv[2], NULL, 0);
	if (errno || nrep < 0 || nrep > 100000) {
		fprintf(stderr, "Error: ``%s'' not in range 0-100000.\n", argv[2]);
		return -2;
	}

	n = argc - 3;
	sizes = calloc(sizeof(int), n);
	if (sizes == NULL) {
		perror("calloc");
		return -3;
	}

	for (i=0; i<n; i++) {
		errno = 0;
		sizes[i] = strtoul(argv[i+3], NULL, 0);
		if (errno || sizes[i] < 0 || sizes[i] > 10000) {
			fprintf(stderr, "Error: ``%s'' not in range 0-10000.\n", argv[i+3]);
			return -4;
		}
	}

	// sort sizes, just in case
	qsort(sizes, n, sizeof(int), cmp_n);

	if (subset(argv[1], nrep, sizes, n) < 0)
		return -5;

	return 0;
}
