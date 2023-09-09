#ifndef __ESTIMATOR__
#define __ESTIMATOR__

#include "lut.h"

#include <stdio.h>
#include <math.h>




void gl_log10(int base, double errate, double *like);

void rescale_likelihood_ratio(double *like);




#endif // __ESTIMATOR__
