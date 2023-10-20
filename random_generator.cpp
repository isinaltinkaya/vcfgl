#include "random_generator.h"
#include "io.h"
#include "shared.h"

#define PI 3.141592654

BetaSampler::BetaSampler(const double mean, const double var, const int seed) {
    double oneOverMean = 1.0 / mean;

    this->alpha = (((1.0 - mean) / var) - (oneOverMean)) * pow(mean, 2);
    this->beta = this->alpha * ((oneOverMean)-1);

    if (alpha <= 0.0) {
        ERROR(
            "Alpha shape parameter of beta distribution is estimated as %f, "
            "which is less than the minimum allowed value of %f. Please use "
            "different --error-rate and --beta-variance values. Current values "
            "are: --error-rate %f --beta-variance %e\n",
            alpha, 0.0, mean, var);
    }
    if (beta <= 0.0) {
        ERROR(
            "Beta shape parameter of beta distribution is estimated as %f, "
            "which is less than the minimum allowed value of %f. Please use "
            "different --error-rate and --beta-variance values. Current values "
            "are: --error-rate %f --beta-variance %e\n",
            beta, 0.0, mean, var);
    }

    fprintf(stderr,
        "\n-> Beta distribution shape parameters are estimated as alpha=%f "
        "and beta=%f (mean=%f, variance=%e)\n",
        this->alpha, this->beta, mean, var);

    fprintf(args->arg_ff, "\n-> Beta distribution shape parameters are estimated as alpha=%f and beta=%f (mean=%f, variance=%e)\n", this->alpha, this->beta, mean, var);

    this->gamma_alpha = new std::gamma_distribution<double>(this->alpha, 1.0);
    this->gamma_beta = new std::gamma_distribution<double>(this->beta, 1.0);

    generator.seed(static_cast<unsigned long>(seed));
}

BetaSampler::~BetaSampler() {
    delete gamma_alpha;
    delete gamma_beta;
}

double BetaSampler::sample() {
    std::gamma_distribution<double> gamma_alpha(alpha, 1.0);
    std::gamma_distribution<double> gamma_beta(beta, 1.0);
    double x = gamma_alpha(generator);
    double y = gamma_beta(generator);

    double ret = (x / (x + y));

    return (ret);
}

int sample_uniform_from_range(int min, int max) {
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
double Poisson(double xm) {
    double lgamma(double xx);
    static double sq, alxm, g, oldm = (-1.0);
    double em, t, y;

    if (xm < 12.0) {
        if (xm != oldm) {
            oldm = xm;
            g = exp(-xm);
        }
        em = -1;
        t = 1.0;
        do {
            ++em;
            t *= drand48();
        } while (t > g);
    } else {
        if (xm != oldm) {
            oldm = xm;
            sq = sqrt(2.0 * xm);
            alxm = log(xm);
            g = xm * alxm - lgamma(xm + 1.0);
        }
        do {
            do {
                y = tan(PI * drand48());
                em = sq * y + xm;
            } while (em < 0.0);
            em = floor(em);
            t = 0.9 * (1.0 + y * y) * exp(em * alxm - lgamma(em + 1.0) - g);
        } while (drand48() > t);
    }
    return em;
}


/// @brief sample_strand - sample a strand
/// @return		0: forward, 1: reverse
int sample_strand() {
    return (drand48() < 0.5 ? 0 : 1);
}

/// @brief sample_tail_distance - sample a tail length for I16 tag
/// @return		int length
int sample_tail_distance() {
    int i;
    if ((i = sample_uniform_from_range(1, 50)) > CAP_DIST) {
        return (CAP_DIST);
    }
    return (i);
}
