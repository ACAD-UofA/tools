/*
 * Copyright (c) 2016,2017 Graham Gower <graham.gower@gmail.com>
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
#include <sys/types.h>
#include <unistd.h>
#include <getopt.h>
#include <errno.h>
#include <ctype.h>

#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>

#define MAX_DP 1024

typedef struct {
	int ignore_transitions;
	char *filter_str;
	char *regions_fn;
	char *regions;
	char *oprefix;
} opt_t;


int
vcfhet(opt_t *opt, char **vcflist, int n)
{
	bcf_srs_t *sr;
	bcf_hdr_t *hdr;
	bcf1_t *rec;
	int ret;
	int i, j;
	int nr;
	int nsamples = 0;
       
	int ndpr_arr = 0;
	int32_t *dpr_arr = NULL; // DPR or AD array

	char ref, alt;
	uint32_t dp_ref, dp_alt;;

	uint64_t het[MAX_DP] = {0,};
	uint64_t nsites[MAX_DP] = {0,};

	sr = bcf_sr_init();
	if (sr == NULL) {
		ret = -1;
		goto err0;
	}

	sr->require_index = 1;
	sr->collapse = COLLAPSE_ANY;

	if (opt->regions_fn && bcf_sr_set_regions(sr, opt->regions_fn, 1) < 0) {
		ret = -2;
		goto err0;
	} else if (opt->regions && bcf_sr_set_regions(sr, opt->regions, 0) < 0) {
		ret = -2;
		goto err0;
	}

	for (i=0; i<n; i++) {
		if (bcf_sr_add_reader(sr, vcflist[i]) == 0) {
			fprintf(stderr, "%s: %s\n", vcflist[i], bcf_sr_strerror(sr->errnum));
			ret = -3;
			goto err0;
		}
		nsamples += bcf_hdr_nsamples(bcf_sr_get_header(sr, i));
	}

	unsigned short xsubi[3] = {31,41,59}; // random state

	while ((nr = bcf_sr_next_line(sr))) {

		ref = '.';
		alt = '.';
		dp_ref = 0;
		dp_alt = 0;

		for (i=0; i<n; i++) {
			hdr = bcf_sr_get_header(sr, i);
			rec = bcf_sr_get_line(sr, i);

			if (rec == NULL || rec->n_allele > 2) {
				// ignore record
				continue;
			}

			if (opt->filter_str && bcf_has_filter(hdr, rec, opt->filter_str) > 0)
				// ignore site
				break;

			//printf("[%s:%d, %d] nr=%d, qual=%lf, n_allele=%d\n", bcf_seqname(hdr, rec), rec->pos+1, i, nr, rec->qual, rec->n_allele);

			if (bcf_get_format_int32(hdr, rec, "DPR", &dpr_arr, &ndpr_arr) < 1) {
				// DPR is deprecated, try AD
				if (bcf_get_format_int32(hdr, rec, "AD", &dpr_arr, &ndpr_arr) < 1) {
					fprintf(stderr, "Error: %s: no FORMAT/DPR nor FORMAT/AD field at %s:%d.\n",
						vcflist[i], bcf_seqname(hdr, rec), rec->pos+1);
					ret = -4;
					goto err1;
				}
			}

			if (ref == '.')
				ref = rec->d.allele[0][0];

			for (j=0; j<bcf_hdr_nsamples(hdr); j++) {
				if (rec->n_allele == 1) {
					dp_ref += dpr_arr[j];
					continue;
				}

				// rec->n_allele == 2
				assert(rec->n_allele == 2);

				// only take SNPs
				if (bcf_get_variant_type(rec, 1) != VCF_SNP)
					break;

				char alt_i = rec->d.allele[1][0];
				if (alt == '.')
					alt = alt_i;
				else if (alt != alt_i) {
					// different ALT genotypes called in different samples
					break;
				}

				dp_ref += dpr_arr[j*2];
				dp_alt += dpr_arr[j*2+1];
			}

			if (j != bcf_hdr_nsamples(hdr))
				// ignore site
				break;
		}

		if (i != n)
			// ignore site
			continue;

		int dp = dp_ref + dp_alt;
		if (dp == 0)
			continue;

		// filter transitions
		if (opt->ignore_transitions &&
		   ((ref == 'C' && alt == 'T') || (alt == 'C' && ref == 'T') ||
		   (ref == 'G' && alt == 'A') || (alt == 'G' && ref == 'A')))
			continue;

		// randomly select two reads
		int gt0 = (nrand48(xsubi) % dp) < dp_ref;
		int gt1 = (nrand48(xsubi) % dp) < dp_ref;

		if (dp < MAX_DP) {
			het[dp] += (gt0 != gt1);
			nsites[dp]++;
		} else {
			fprintf(stderr, "Warning: dp(%d) >= MAX_DP(%d) at %s:%d, consider increasing MAX_DP\n",
					dp, MAX_DP, bcf_seqname(hdr, rec), rec->pos+1);
		}
	}

	printf("DP\tNsites\tHets\n");
	for (i=1; i<MAX_DP; i++) {
		printf("%d\t%jd\t%jd\n", i, (uintmax_t)nsites[i], (uintmax_t)het[i]);
		if (nsites[i]==0)
			break;
	}

	ret = 0;
err1:
	if (dpr_arr != NULL)
		free(dpr_arr);

        bcf_sr_destroy(sr);
err0:
	return ret;
}

void
usage(char *argv0)
{
	fprintf(stderr, "usage: %s [...] file1.vcf [... fileN.vcf]\n", argv0);
	fprintf(stderr, "   -r STR           Comma separated list of regions to include []\n");
	fprintf(stderr, "   -R FILE          Bed file listing regions to include []\n");
	fprintf(stderr, "   -t               Exclude transitions [no]\n");
	fprintf(stderr, "   -F STR           Ignore sites with STR in any file's FILTER column [SnpGap]\n");
	//fprintf(stderr, "   -o STR           Output file prefix [out]\n");
	exit(1);
}

int
main(int argc, char **argv)
{
	opt_t opt;
	int c;

	opt.filter_str = NULL;
	//opt.filter_str = "SnpGap";
	opt.regions_fn = NULL;
	opt.regions = NULL;
	opt.ignore_transitions = 0;
	opt.oprefix = "out";

	while ((c = getopt(argc, argv, "r:R:o:F:t")) != -1) {
		switch (c) {
			case 't':
				opt.ignore_transitions = 1;
				break;
			case 'R':
				opt.regions_fn = optarg;
				break;
			case 'r':
				opt.regions = optarg;
				break;
			case 'o':
				opt.oprefix = optarg;
				break;
			case 'F':
				opt.filter_str = optarg;
				break;
			default:
				usage(argv[0]);
		}
	}

	if (argc-optind < 1) {
		usage(argv[0]);
	}

	if (opt.regions && opt.regions_fn) {
		fprintf(stderr, "-R and -r are mutually exclusive.");
		usage(argv[0]);
	}

	if (vcfhet(&opt, argv+optind, argc-optind) < 0)
		return -1;

	return 0;
}

