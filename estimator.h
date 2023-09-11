#ifndef __ESTIMATOR__
#define __ESTIMATOR__

#include <math.h>


#include "shared.h"




void gl_log10(int base, double errate, double *like);
void gl_ln(int base, double errate, double *like);

void rescale_likelihood_ratio(double *like);
void rescale_likelihood_ratio(float *like);

void rescale_likelihood_ratio(float* like, const int size);



#endif // __ESTIMATOR__
