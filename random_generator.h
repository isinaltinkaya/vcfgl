#ifndef __RANDOM_GENERATOR__
#define __RANDOM_GENERATOR__

#include <random>
#include "shared.h"

#include <stdio.h>

#include <inttypes.h>
#include <math.h>

#include <stdlib.h>

struct BetaSampler
{

    double alpha, beta;
    std::mt19937 generator;
    std::gamma_distribution<double> *gamma_alpha;
    std::gamma_distribution<double> *gamma_beta;

    BetaSampler(const double mean, const double var, const int seed);
    ~BetaSampler();

    double sample();
};

int sample_uniform_from_range(int min, int max);

double Poisson(double xm);

/// @brief sample_strand - sample a strand
/// @return		0: forward, 1: reverse
inline int sample_strand()
{
	return (drand48() < 0.5 ? 0 : 1);
}

/// @brief sample_tail_distance - sample a tail length for I16 tag
/// @return		int length
inline int sample_tail_distance()
{
	int i;
	if ((i = sample_uniform_from_range(1, 50)) > CAP_DIST)
	{
		return (CAP_DIST);
	}
	return (i);
}

#endif
