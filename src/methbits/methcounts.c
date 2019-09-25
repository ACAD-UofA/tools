#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <math.h>

#define COVMAX 1024
#define LOG2COVMAX 10

int
parse_pileometh(char *fn)
{
	int ret;
	FILE *fp;
	char *line, *p;
	size_t buflen;
	ssize_t n_bytes;
	uint64_t lineno;
	unsigned int i, m, u;

	uint64_t *counts;

	fp = fopen(fn, "r");
	if (fp == NULL) {
		fprintf(stderr, "fopen: %s: %s\n", fn, strerror(errno));
		ret = -1;
		goto err0;
	}

	counts = calloc(COVMAX*COVMAX, sizeof(uint64_t));
	if (counts == NULL) {
		perror("calloc");
		ret = -2;
		goto err1;
	}

	line = NULL;
	buflen = 0;

	// get header, ignore
	if ((n_bytes = getline(&line, &buflen, fp)) == -1) {
		fprintf(stderr, "getline: %s:1: %s\n", fn, strerror(errno));
		ret = -3;
		goto err2;
	}
	if (strcmp("track type", line))
		// that wasn't a header...
		rewind(fp);

	errno = 0;

	for (lineno=2; ; lineno++) {
		if ((n_bytes = getline(&line, &buflen, fp)) == -1) {
			if (feof(fp))
				break;
			fprintf(stderr, "getline: %s:%jd: %s\n", fn, lineno,
					strerror(errno));
			ret = -3;
			goto err3;
		}

		if (n_bytes == 0)
			// eof
			break;

		p = line;

		if (p[n_bytes-1] == '\n')
			p[n_bytes-1] = '\0';

// skip to next column
#define next(x) \
		while (*x && *x != ' ' && *x != '\t') x++; \
		while (*x && (*x == ' ' || *x == '\t')) x++;

		// columns are: chr start end meth% methylated unmethylated

		next(p);
		if (!*p)
			// empty line
			continue;

		next(p);
		next(p);
		next(p);

		m = strtoul(p, &p, 0);
		if (errno) {
			fprintf(stderr, "strtoul:1: %s:%jd: %s\n", fn, lineno,
					strerror(errno));
			ret = -4;
			goto err3;
		}

		next(p);
		u = strtoul(p, &p, 0);
		if (errno) {
			fprintf(stderr, "strtoul:2: %s:%jd: %s\n", fn, lineno,
					strerror(errno));
			ret = -5;
			goto err3;
		}

		if (m+u >= COVMAX)
			continue;

		counts[m<<LOG2COVMAX | u]++;
	}

	printf("methylated\tunmethylated\tcount\n");
	for (i=0; i<COVMAX*COVMAX; i++) {
		if (counts[i] == 0)
			continue;
		m = i >> LOG2COVMAX;
		u = i & (COVMAX-1);
		printf("%d\t%d\t%jd\n", m, u, counts[i]);
	}

	ret = 0;
err3:
	if (buflen)
		free(line);
err2:
	free(counts);
err1:
	fclose(fp);
err0:
	return ret;
}

int
parse_methylkit(char *fn)
{
	int ret;
	FILE *fp;
	char *line, *p;
	size_t buflen;
	ssize_t n_bytes;
	uint64_t lineno;
	unsigned int i, coverage, m, u;
	double freqC;

	uint64_t *counts;

	fp = fopen(fn, "r");
	if (fp == NULL) {
		fprintf(stderr, "fopen: %s: %s\n", fn, strerror(errno));
		ret = -1;
		goto err0;
	}

	counts = calloc(COVMAX*COVMAX, sizeof(uint64_t));
	if (counts == NULL) {
		perror("calloc");
		ret = -2;
		goto err1;
	}

	line = NULL;
	buflen = 0;

	// get header, ignore
	if ((n_bytes = getline(&line, &buflen, fp)) == -1) {
		fprintf(stderr, "getline: %s:1: %s\n", fn, strerror(errno));
		ret = -3;
		goto err2;
	}

	errno = 0;

	for (lineno=2; ; lineno++) {
		if ((n_bytes = getline(&line, &buflen, fp)) == -1) {
			if (feof(fp))
				break;
			fprintf(stderr, "getline: %s:%jd: %s\n", fn, lineno,
					strerror(errno));
			ret = -3;
			goto err3;
		}

		if (n_bytes == 0)
			// eof
			break;

		p = line;

		if (p[n_bytes-1] == '\n')
			p[n_bytes-1] = '\0';

// skip to next column
#define next(x) \
		while (*x && *x != ' ' && *x != '\t') x++; \
		while (*x && (*x == ' ' || *x == '\t')) x++;

		// columns are: chrBase chr base strand coverage freqC freqT

		next(p);
		if (!*p)
			// empty line
			continue;

		next(p);
		next(p);
		next(p);

		coverage = strtoul(p, &p, 0);
		if (errno) {
			fprintf(stderr, "strtoul: %s:%jd: %s\n", fn, lineno,
					strerror(errno));
			ret = -4;
			goto err3;
		}

		if (coverage >= COVMAX)
			continue;

		next(p);
		freqC = strtod(p, &p);
		if (errno) {
			fprintf(stderr, "strtod: %s:%jd: %s\n", fn, lineno,
					strerror(errno));
			ret = -5;
			goto err3;
		}

		m = (unsigned int) round(0.01*coverage*freqC);
		u = coverage -m;

		counts[m<<LOG2COVMAX | u]++;
	}

	printf("methylated\tunmethylated\tcount\n");
	for (i=0; i<COVMAX*COVMAX; i++) {
		if (counts[i] == 0)
			continue;
		m = i >> LOG2COVMAX;
		u = i & (COVMAX-1);
		printf("%d\t%d\t%jd\n", m, u, counts[i]);
	}

	ret = 0;
err3:
	if (buflen)
		free(line);
err2:
	free(counts);
err1:
	fclose(fp);
err0:
	return ret;
}

int
main(int argc, char **argv)
{
	return parse_pileometh("/dev/stdin");
	//return parse_methylkit("/dev/stdin");
}

