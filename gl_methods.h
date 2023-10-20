#ifndef __GL_METHODS__
#define __GL_METHODS__

#include <math.h>

#include "shared.h"

void calc_gl_log10_method1(int base, double e, double* acgt_gls);
void calc_gl_log10_method1_precalc(int base, double e, double* acgt_gls);

void gl_log10(int base, double error_rate, double* like);
void gl_ln(int base, double error_rate, double* like);

void rescale_likelihood_ratio(double* like);
void rescale_likelihood_ratio(float* like);

void rescale_likelihood_ratio(float* like, const int size);

#endif  // __GL_METHODS__
