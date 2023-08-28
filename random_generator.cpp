#include "random_generator.h"

#include <stdio.h>

#include <inttypes.h>
#include <math.h>

#define PI 3.141592654

int sample_uniform_from_range(int min, int max){
   return (min + rand() / (RAND_MAX / (max - min + 1) + 1));
}





// Poisson random deviates
//
// From tsk msToGlf.cpp
// See NRC++ e2 p.298
//
// returns em= a random deviate drawn from a poisson dist
// of mean xm, using drand48() as a source of uniform
// random deviates with a srand48 pseudo-random number
// initializer
double Poisson(double xm)
{
	double lgamma(double xx);
	static double sq,alxm,g,oldm=(-1.0);
	double em, t, y;

	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em=-1;
		t=1.0;
		do {
			++em;
			t *=drand48();
		} while (t>g);
	}
	else {
		if (xm!=oldm) {
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-lgamma(xm+1.0);
		}
		do {
			do {
				y=tan(PI*drand48());
				em=sq*y+xm;
			} while (em< 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-lgamma(em+1.0)-g);
		} while (drand48()>t);
	}
	return em;
}


