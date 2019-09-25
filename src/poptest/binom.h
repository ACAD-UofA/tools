/*
 * Copyright (c) 2015 Graham Gower <graham.gower@gmail.com>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 */

#ifndef BINOM_H
#define BINOM_H

double lchoose(double k, double n);
double binom_lpmf(double k, double n, double p);
double int_binom_wrt_p(int k, int n);

#endif
