#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include "svgfig.h"

typedef struct {
	char *fill;
	char *stroke;
	float stroke_width;
} style_t;

svgfig_t *
svgfig_new(char *fn, int w, int h)
{
	svgfig_t *f;
	f = malloc(sizeof(svgfig_t));
	if (f == NULL) {
		perror("malloc");
		return NULL;
	}
	f->fn = fn;
	f->w = w;
	f->h = h;

	f->fp = fopen(fn, "w");
	if (f->fp == NULL) {
		fprintf(stderr, "%s: %s\n", fn, strerror(errno));
	}

	fprintf(f->fp, "<svg"
			" viewBox='0 0 1 1' preserveAspectRatio='none'"
			" width='%d' height='%d'"
			" xmlns='http://www.w3.org/2000/svg'"
			" xmlns:xlink='http://www.w3.org/1999/xlink'"
			">\n", w, h);

	return f;
}

void
svgfig_close(svgfig_t *f)
{
	fprintf(f->fp, "</svg>\n");
	fclose(f->fp);
	free(f);
}

// styling via presentation attributes
void
svg_style(svgfig_t *f, style_t *style)
{
	if (!style)
		return;

	fprintf(f->fp, " fill='%s'", style->fill ? style->fill : "none");
	if (style->stroke && style->stroke_width > 0) {
		fprintf(f->fp, " stroke-width='%g' stroke='%s'",
			style->stroke_width, style->stroke);
	}
}

void
svg_rect(svgfig_t *f, float x, float y, float w, float h, style_t *style)
{
	fprintf(f->fp, "<rect x='%g' y='%g' width='%g' height='%g'",
		x, y, w, h);
	svg_style(f, style);
	fprintf(f->fp, " vector-effect='non-scaling-stroke'");
	fprintf(f->fp, " />\n");
}

void
svg_line(svgfig_t *f, float x1, float y1, float x2, float y2, style_t *style)
{
	fprintf(f->fp, "<line x1='%g' y1='%g' x2='%g' y2='%g'",
			x1, y1, x2, y2);
	svg_style(f, style);
	fprintf(f->fp, " vector-effect='non-scaling-stroke'");
	fprintf(f->fp, " />\n");
}

void
svg_text(svgfig_t *f, char *text, float x, float y, style_t *style)
{
	fprintf(f->fp, "<text x='%g' y='%g'", x, y);
	svg_style(f, style);
	//fprintf(f->fp, " vector-effect='non-scaling-stroke'");
	fprintf(f->fp, " font-size='0.1%%' text-anchor='middle' dominant-baseline='central' font-family='Arial, sans-serif'");
	fprintf(f->fp, ">%s</text>\n", text);
}

void
extent(double *x, size_t n, double *min, double *max)
{
	int i;
	*min = DBL_MAX;
	*max = DBL_MIN;
	for (i=0; i<n; i++) {
		if (x[i] < *min)
			*min = x[i];
		if (x[i] > *max)
			*max = x[i];
	}
}

void
svgfig_bar(svgfig_t *f, double *x, double *y, size_t n, double width)
{
	int i;
	double w, h;
	double xlim[2], ylim[2];

	/*
	 * bounding box for the plot
	 */
	fprintf(f->fp, "<g transform='translate(0.1 0.1) scale(0.8 0.8)'>\n");
	svg_rect(f, 0, 0, 1, 1, &(style_t){.stroke="black",.stroke_width=1});

	extent(x,n,xlim,xlim+1);
	extent(y,n,ylim,ylim+1);

	w = xlim[1]-xlim[0] + width;
	h = ylim[1];

	/*
	 * Transform and scale the plotting region for cartesian coords.
	 * This works in Firefox, but remains upside down in Chrome.
	 */
	fprintf(f->fp, "<svg viewBox='%g %g %g %g'"
			" transform='translate(0 1) scale(1 -1)'"
			" preserveAspectRatio='none'",
			xlim[0], 0.0, w, h);
	svg_style(f, &(style_t){.fill="blue",.stroke="black",.stroke_width=1});
	fprintf(f->fp, ">\n");

	/*
	 * draw the bars
	 */
	for (i=0; i<n; i++) {
		svg_rect(f, xlim[0]+x[i], 0, width, y[i], NULL);
	}
	fprintf(f->fp, "</svg>\n</g>\n");

	/*
	 * tick marks and labels
	 */
	char buf[128];
	double xtick, ytick;

	fprintf(f->fp, "<g>\n");
	for (i=0; i<n; i++) {
		xtick = 0.1 + 0.8*(xlim[0]+x[i])/w;
		svg_line(f, xtick, 0.9, xtick, 0.92, &(style_t){.stroke="black",.stroke_width=1});
		snprintf(buf, 127, "%g", x[i]);
		svg_text(f, buf, xtick, 0.95, &(style_t){.fill="black"});

		ytick = 0.9 -0.8*(double)i/n;
		svg_line(f, 0.09, ytick, 0.1, ytick, &(style_t){.stroke="black",.stroke_width=1});
		snprintf(buf, 127, "%g", i*h/n);
		svg_text(f, buf, 0.05, ytick, &(style_t){.fill="black"});
	}
	xtick = 0.1 + 0.8*(xlim[0]+x[i-1]+width)/w;
	svg_line(f, xtick, 0.9, xtick, 0.92, &(style_t){.stroke="black",.stroke_width=1});
	snprintf(buf, 127, "%g", x[i-1]+width);
	svg_text(f, buf, xtick, 0.95, &(style_t){.fill="black"});

	ytick = 0.1;
	svg_line(f, 0.09, ytick, 0.1, ytick, &(style_t){.stroke="black",.stroke_width=1});
	snprintf(buf, 127, "%g", h);
	svg_text(f, buf, 0.05, ytick, &(style_t){.fill="black"});
	fprintf(f->fp, "</g>\n");
}

void
svgfig_hist(svgfig_t *f, double *x, size_t n, int nbins)
{
	int i;
	double min, max;
	double range;
	double *bins;
	double *x2;

	if (nbins == 0)
		return;

	bins = calloc(nbins, sizeof(double));
	x2 = calloc(nbins, sizeof(double));
	if (bins == NULL || x2 == NULL) {
		perror("calloc");
		return;
	}

	extent(x, n, &min, &max);
	range = (max - min) / nbins;

	for (i=0; i<n; i++) {
		int bin = (x[i] - min) / range;
		bins[bin]++;
	}

	for (i=0; i<nbins; i++) {
		x2[i] = min + i*range;
	}

	svgfig_bar(f, x2, bins, nbins, range);
	free(bins);
	free(x2);
}
