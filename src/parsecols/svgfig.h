#ifndef __SVGFIG_H

typedef struct {
	char *fn;
	FILE *fp;
	int w, h;
} svgfig_t;

svgfig_t *svgfig_new(char *, int, int);
void svgfig_close(svgfig_t *);
void svgfig_hist(svgfig_t *, double *, size_t, int);

#endif // __SVGFIG_H
