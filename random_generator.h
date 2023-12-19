#ifndef __RANDOM_GENERATOR__
#define __RANDOM_GENERATOR__

#include <random> // std::mt19937, std::gamma_distribution

extern unsigned short int rng1_seeder[3];
extern unsigned short int rng1_seeder_save[3];

struct BetaSampler {
    double alpha, beta;

    std::mt19937 generator;

    std::gamma_distribution<double>* gamma_alpha;
    std::gamma_distribution<double>* gamma_beta;

    BetaSampler(const double mean, const double var, const int seed);
    ~BetaSampler();

    double sample();
};


double sample_Poisson_rng0(double xm);
double sample_Poisson_rng1(double xm);
double sample_Poisson_rng1_rewind(double xm);

double sample_uniform_rng0(void);
double sample_uniform_rng1(void);
double sample_uniform_rng2(void);

inline int sample_from_range_rng2(int min, int max) {
    return (min + rand() / (RAND_MAX / (max - min + 1) + 1));
}

inline int sample_from_range_rng0(int min, int max) {
    return (min + sample_uniform_rng0() * (max - min + 1));
}

#endif
