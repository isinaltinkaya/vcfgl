#include "random_generator.h"
#include "io.h"

NormalSampler* NormalSampler_init(const double mean, const double var, const int seed) {

    NormalSampler* norm = (NormalSampler*)malloc(sizeof(NormalSampler));

    norm->mu = mean;
    norm->sigma = var;

    return(norm);
}

void NormalSampler_destroy(NormalSampler* norm) {
    free(norm);
    norm = NULL;
}


Gamma1Sampler* Gamma1Sampler_init(const double shape) {

    DEVASSERT(shape > 0.0);

    Gamma1Sampler* gamma = (Gamma1Sampler*)malloc(sizeof(Gamma1Sampler));

    gamma->alpha = shape;
    gamma->alpha_changed = false;
    if (gamma->alpha < 1.0) {
        gamma->alpha += 1.0;
        gamma->alpha_changed = true;
    }

    gamma->a1 = gamma->alpha - 1.0 / 3.0;
    gamma->a2 = 1.0 / sqrt(9. * gamma->a1);

    //TODO write simpler normalsampler for this specific instance
    gamma->normalSampler = NormalSampler_init(0.0, 1.0, 0);

    return(gamma);
}

void Gamma1Sampler_destroy(Gamma1Sampler* gamma) {
    NormalSampler_destroy(gamma->normalSampler);
    free(gamma);
    gamma = NULL;
}

double Gamma1Sampler::sample(void) {
    double u, v, x, xsq;
    do {
        do {
            x = this->normalSampler->sample();
            v = 1.0 + this->a2 * x;
        } while (v <= 0.0);
        v = v * v * v;
        u = sample_uniform_rng2();
        xsq = x * x;
    } while (u > 1.0 - 0.0331 * (xsq * xsq) &&
        log(u) > 0.5 * xsq + this->a1 * (1.0 - v + log(v)));
    if (this->alpha_changed) {
        while ((u = sample_uniform_rng2()) == 0.0);
        return (pow(u, 1.0 / this->old_alpha) * this->a1 * v);
    } else {
        return (this->a1 * v);
    }

}

GammaSampler* GammaSampler_init(const double shape, const double rate) {

    DEVASSERT(shape > 0.0);
    DEVASSERT(rate > 0.0);

    GammaSampler* gamma = (GammaSampler*)malloc(sizeof(GammaSampler));

    gamma->alpha = shape;
    gamma->alpha_changed = false;
    if (gamma->alpha < 1.0) {
        gamma->alpha += 1.0;
        gamma->alpha_changed = true;
    }

    gamma->beta = rate;
    gamma->a1 = gamma->alpha - 1.0 / 3.0;
    gamma->a2 = 1.0 / sqrt(9. * gamma->a1);
    gamma->normalSampler = NormalSampler_init(0.0, 1.0, 0);

    return(gamma);
}

void GammaSampler_destroy(GammaSampler* gamma) {
    NormalSampler_destroy(gamma->normalSampler);
    free(gamma);
    gamma = NULL;
}


double GammaSampler::sample(void) {
    double u, v, x, xsq;
    do {
        do {
            x = this->normalSampler->sample();
            v = 1.0 + this->a2 * x;
        } while (v <= 0.0);
        v = v * v * v;
        u = sample_uniform_rng2();
        xsq = x * x;
    } while (u > 1.0 - 0.0331 * (xsq * xsq) &&
        log(u) > 0.5 * xsq + this->a1 * (1.0 - v + log(v)));
    if (this->alpha_changed) {
        while ((u = sample_uniform_rng2()) == 0.0);
        return (pow(u, 1.0 / this->old_alpha) * this->a1 * v / this->beta);
    } else {
        return (this->a1 * v / this->beta);
    }

}


// todo set ptr once instead of this check
double PoissonSampler::sample(void) {
    if (this->st12) {
        return(sample_st12());
    } else {
        return(sample_bteq12());
    }
}

double PoissonSampler::sample_st12(void) {
    double em = -1.0;
    double t = 1.0;
    do {
        ++em;
        t *= erand48(rng1_seeder);
    } while (t > (this->g));
    return(em);
}

double PoissonSampler::sample_bteq12(void) {
    double em, t, y;
    do {
        do {
            y = tan(PI * erand48(rng1_seeder));
            em = (this->sq) * y + (this->lm);
        } while (em < 0.0);
        em = floor(em);
        t = 0.9 * (1.0 + y * y) * exp(em * (this->alxm) - gamma_ln(em + 1.0) - (this->g));
    } while (erand48(rng1_seeder) > t);
    return(em);
}

PoissonSampler* PoissonSampler_init(const double lambda) {

    PoissonSampler* pois = (PoissonSampler*)malloc(sizeof(PoissonSampler));

    pois->lm = lambda;

    // init
    pois->sq = -1.0;
    pois->alxm = -1.0;
    pois->g = -1.0;
    pois->st12 = true;

    if (lambda < 12.0) {
        pois->g = exp(-(pois->lm));
    } else {
        pois->st12 = false;
        pois->sq = sqrt(2.0 * pois->lm);
        pois->alxm = log(pois->lm);
        pois->g = pois->lm * pois->alxm - gamma_ln(pois->lm + 1.0);
    }
    return(pois);
}



#if __USE_STD_BETA__ == 1
BetaSampler::BetaSampler(const double mean, const double var, const int seed) {

    ASSERT(mean > 0.0);
    ASSERT(var > 0.0);
    ASSERT(mean < 1.0);

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

    fprintf(args->arg_fp, "\n-> Beta distribution shape parameters are estimated as alpha=%f and beta=%f (mean=%f, variance=%e)\n", this->alpha, this->beta, mean, var);

    this->gamma_x = new std::gamma_distribution<double>(this->alpha, 1.0);
    this->gamma_y = new std::gamma_distribution<double>(this->beta, 1.0);

    generator.seed(seed);
}

BetaSampler::~BetaSampler() {
    delete gamma_x;
    delete gamma_y;
}

double BetaSampler::sample(void) {
    std::gamma_distribution<double> gamma_x(alpha, 1.0);
    std::gamma_distribution<double> gamma_y(beta, 1.0);
    double x = gamma_x(generator);
    double y = gamma_y(generator);

    double ret = (x / (x + y));
    ASSERT(ret >= 0.0);
    ASSERT(ret <= 1.0);

    return (ret);
}

#else

BetaSampler* BetaSampler_init(const double mean, const double var, const int seed) {

    ASSERT(mean > 0.0);
    ASSERT(var > 0.0);
    ASSERT(mean < 1.0);

    double oneOverMean = 1.0 / mean;

    BetaSampler* beta = (BetaSampler*)malloc(sizeof(BetaSampler));

    beta->alpha = (((1.0 - mean) / var) - (oneOverMean)) * pow(mean, 2);
    beta->beta = beta->alpha * ((oneOverMean)-1);

    if (beta->alpha <= 0.0) {
        ERROR(
            "Alpha shape parameter of beta distribution is estimated as %f, "
            "which is less than the minimum allowed value of %f. Please use "
            "different --error-rate and --beta-variance values. Current values "
            "are: --error-rate %f --beta-variance %e\n",
            beta->alpha, 0.0, mean, var);
    }
    if (beta->beta <= 0.0) {
        ERROR(
            "Beta shape parameter of beta distribution is estimated as %f, "
            "which is less than the minimum allowed value of %f. Please use "
            "different --error-rate and --beta-variance values. Current values "
            "are: --error-rate %f --beta-variance %e\n",
            beta->beta, 0.0, mean, var);
    }

    fprintf(stderr,
        "\n-> Beta distribution shape parameters are estimated as alpha=%f "
        "and beta=%f (mean=%f, variance=%e)\n",
        beta->alpha, beta->beta, mean, var);

    fprintf(args->arg_fp, "\n-> Beta distribution shape parameters are estimated as alpha=%f and beta=%f (mean=%f, variance=%e)\n", beta->alpha, beta->beta, mean, var);


    // beta->gamma_x = GammaSampler_init(beta->alpha, 1.0);
    // beta->gamma_y = GammaSampler_init(beta->beta, 1.0);
    beta->gamma_x = Gamma1Sampler_init(beta->alpha);
    beta->gamma_y = Gamma1Sampler_init(beta->beta);

    return(beta);
}

void BetaSampler_destroy(BetaSampler* beta) {
    Gamma1Sampler_destroy(beta->gamma_x);
    Gamma1Sampler_destroy(beta->gamma_y);
    free(beta);
    beta = NULL;
}

double BetaSampler::sample(void) {

    double x = this->gamma_x->sample();
    double y = this->gamma_y->sample();



    double ret = (x / (x + y));

    DEVASSERT(ret >= 0.0);
    DEVASSERT(ret <= 1.0);

    return (ret);
}

#endif