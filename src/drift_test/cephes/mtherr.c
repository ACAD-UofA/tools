#include "mconf.h"

int merror = 0;

static char *c2str[] = {
	[0] = "unknown error",
	[DOMAIN] = "argument domain error",
	[SING] = "argument singularity",
	[OVERFLOW] = "overflow range error",
	[UNDERFLOW] = "underflow range error",
	[TLOSS] = "total loss of precision",
	[PLOSS] = "partial loss of precision",
	[TOOMANY] = "too many iterations",
};

int
mtherr(char *name, int code)
{
	if (code < 0 || code > 8)
		code = 0;
	merror = code;
	if (code != PLOSS)
		fprintf(stderr, "%s: %s\n", name, c2str[code]);
	return 0;
}
