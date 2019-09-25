#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <stdarg.h>

#include "svgfig.h"

#define max(a,b) ((a)>(b)?(a):(b))

typedef struct {
	void **x;
	int col;
	enum {
		INT=1,
		DOUBLE,
		STR,
	} type;
} type_t;

#define NULLTYPE ((type_t){0,})

int
t_compar(const void *v1, const void *v2)
{
	const type_t *t1 = v1;
	const type_t *t2 = v2;
	return t1->col - t2->col;
}

int
parsecols(char *fn, int skip, uint64_t *nlines, ...)
{
	int ret;
	ssize_t n_bytes;
	uint64_t lineno = 0;
	uint64_t memlines = 0;
	FILE *fp;
	char *buf = NULL;
	size_t buflen = 0;

	int col, colmax = -1;
	char *p;

	int i, j;
	type_t *typelist = NULL, *t;
	size_t n_types = 0;
	va_list va;

	*nlines = 0;

	va_start(va, nlines);
	for (;;) {
		void *tmp;
		type_t tt;

		tt = va_arg(va, type_t);
		if (!memcmp(&tt, &NULLTYPE, sizeof(type_t)))
			break;

		if (tt.col == colmax) {
			fprintf(stderr, "Duplicate type specifications for column %d\n", colmax);
			ret = -2;
			goto err0;
		}

		colmax = max(tt.col, colmax);
		(*tt.x) = NULL;

		n_types++;
		tmp = realloc(typelist, sizeof(type_t)*n_types);
		if (tmp == NULL) {
			perror("realloc");
			ret = -1;
			goto err0;
		}
		typelist = tmp;
		typelist[n_types-1] = tt;
	}
	va_end(va);

	if (n_types == 0)
		return 0;

	qsort(typelist, n_types, sizeof(type_t), t_compar);

	fp = fopen(fn, "r");
	if (fp == NULL) {
		fprintf(stderr, "%s: %s\n", fn, strerror(errno));
		ret = -3;
		goto err0;
	}

	j = 0;
	for (;;) {
		lineno++;

		n_bytes = getline(&buf, &buflen, fp);
		if (n_bytes == -1) {
			if (feof(fp))
				break;
			fprintf(stderr, "%s:%jd: %s\n",
					fn, lineno, strerror(errno));
			ret = -4;
			goto err1;
		}

		if (skip) {
			skip--;
			continue;
		}

		if (buf[n_bytes-1] == '\n')
			buf[n_bytes-1] = 0;

		if (lineno-skip > memlines) {
			void *tmp;
			memlines += 1024;
			for (i=0; i<n_types; i++) {
				size_t size;
				t = typelist +i;
				switch (t->type) {
					case INT:
						size = sizeof(int);
						break;
					case DOUBLE:
						size = sizeof(double);
						break;
					case STR:
						size = sizeof(char*);
						break;
					default:
						fprintf(stderr, "unknown type\n");
						abort();
				}
				tmp = realloc(*(t->x), memlines*size);
				if (tmp == NULL) {
					free(*(t->x));
					fprintf(stderr, "%s:%jd: %s\n", fn, lineno, strerror(errno));
					ret = -5;
					goto err1;
				}
				*(t->x) = tmp;
			}
		}

		p = buf;
		t = typelist;
		for (col=0; col<=colmax; col++) {
			char *s = p;
			while (*p != 0 && *p != '\t' && *p != ' ')
				p++;
			if (p == s)
				break;

			while (*p != 0 && (*p == '\t' || * p== ' '))
				*(p++) = 0;

			if (col < t->col)
				continue;

			switch (t->type) {
				case INT:
					{
					int *x = *(int**)t->x;
					x[j] = atoi(s);
					}
					break;
				case DOUBLE:
					{
					double *x = *(double**)t->x;
					x[j] = atof(s);
					}
					break;
				case STR:
					{
					char **x = *(char***)t->x;
					x[j] = strdup(s);
					}
					break;
			}
			t++;
		}

		j++;

		if (col <= colmax) {
			fprintf(stderr, "%s:%jd: need >=%d columns, got %d\n",
					fn, lineno, colmax, col);
			ret = -4;
			goto err1;
		}
	}


	*nlines = j;
	ret = 0;
err1:
	fclose(fp);
err0:
	free(typelist);
	return ret;
}


int
main(int argc, char **argv)
{
	int *pos, *m, *u;
	double *mprop;
	uint64_t n;
	int i;

	if (argc != 2) {
		fprintf(stderr, "usage: %s file.txt\n", argv[0]);
		return -1;
	}

	parsecols(argv[1], 0, &n,
			(type_t){(void**)&pos, 1, INT},
			(type_t){(void**)&m, 4, INT},
			(type_t){(void**)&u, 5, INT},
			NULLTYPE
			);

	mprop = malloc(sizeof(double)*n);
	if (mprop == NULL) {
		perror("malloc");
		return -2;
	}

	for (i=0; i<n; i++) {
		mprop[i] = (double)m[i]/(m[i]+u[i]);
		//printf("%d, %lf\n", pos[i], mprop[i]);
	}

	//int ratio[2] = {4,3};
	int ratio[2] = {16,9};

	svgfig_t *f = svgfig_new("a.svg", ratio[0]*80, ratio[1]*80);
	svgfig_hist(f, mprop, n, 10);
	svgfig_close(f);

	return 0;
}
