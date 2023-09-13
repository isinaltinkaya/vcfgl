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

#endif
