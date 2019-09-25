/*
 * Print invervals of CpG Islands.
 * See: Gardiner-Garden and Frommer, 1987, doi:10.1016/0022-2836(87)90689-9
 *
 * Copyright (c) 2017 Graham Gower <graham.gower@gmail.com>
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
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <fcntl.h>
#include <ctype.h>
#include <errno.h>

#include "kseq.h"
KSEQ_INIT(int, read);

#define WIN 100
#define MINSIZE 200
#define OE_THRES (0.6/WIN)

int
cgi2bed(char *fn)
{
	kseq_t *ks;
	int fd;
	int ret;
	int len;

	uint64_t i, start;
	unsigned int c, g, cpg;
	int island;

	fd = open(fn, O_RDONLY);
	if (fd == -1) {
		fprintf(stderr, "%s: %s\n", fn, strerror(errno));
		ret = -1;
		goto err0;
	}

	ks = kseq_init(fd);

	while ((len = kseq_read(ks)) >= 0) {
		const char *s = ks->seq.s;
		char s0, s1;
		char e0, e1;

		if (ks->seq.l < MINSIZE)
			continue;

		c = g = cpg = 0;

		for (i=0; i<WIN-1; i++) {
			s0 = toupper(s[i]);
			s1 = toupper(s[i+1]);
			switch (s0) {
			case 'C':
				c++;
				if (s1 == 'G')
					cpg++;
				break;
			case 'G':
				g++;
				break;
			}
		}

		if (s1 == 'C')
			c++;
		else if (s1 == 'G')
			g++;

		if (c+g > WIN/2 && cpg && ((double)cpg/(c*g) > OE_THRES)) {
			island = 1;
			start = 0;
		} else
			island = 0;

		for (i=0; i<ks->seq.l-WIN; i++) {
			s0 = toupper(s[i]);
			s1 = toupper(s[i+1]);
			e0 = toupper(s[i+WIN-1]);
			e1 = toupper(s[i+WIN]);

			switch (s0) {
			case 'C':
				c--;
				if (s1 == 'G')
					cpg--;
				break;
			case 'G':
				g--;
				break;
			}

			switch (e1) {
			case 'C':
				c++;
				break;
			case 'G':
				g++;
				if (e0 == 'C')
					cpg++;
				break;
			}

			if (c+g > WIN/2 && cpg && ((double)cpg/(c*g) > OE_THRES)) {
				if (!island) {
					island = 1;
					start = i;
				}
			} else {
				if (island) {
					island = 0;
					if (i+WIN-start >= MINSIZE)
						printf("%s\t%jd\t%jd\n", ks->name.s, start+1, i+WIN);
				}
			}
		}

		if (island && i+WIN-start >= MINSIZE)
			printf("%s\t%jd\t%jd\n", ks->name.s, start+1, i+WIN);
	}

	if (len != -1) {
		fprintf(stderr, "%s: unexpected end of file.\n", fn);
	}

	ret = 0;

	kseq_destroy(ks);
	if (fd != 0)
		close(fd);
err0:
	return ret;
}


int
main(int argc, char **argv)
{
	char *fn;

	if (argc == 1)
		fn = "/dev/stdin";
	else if (argc == 2)
		fn = argv[1];
	else {
		fprintf(stderr, "usage: %s ref.fa\n", argv[0]);
		return -1;
	}

	return cgi2bed(fn);
}
