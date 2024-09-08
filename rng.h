#ifndef __RANDOM_GENERATOR__
#define __RANDOM_GENERATOR__

#include <random> // std::mt19937, std::gamma_distribution

#include "shared.h"

#define sample_uniform_rng0() drand48()
#define sample_uniform_rng1() erand48(rng1_seeder)
#define sample_uniform_rng2() erand48(rng2_seeder)
#define sample_uniform_rng_rand ((double)rand() / (double)RAND_MAX + 1.0)
#define sample_from_range_rng_rand(min,max) ( (min) + (rand()) / (RAND_MAX / ((max) - (min) + 1) +1));

extern unsigned short int rng1_seeder[3];
extern unsigned short int rng1_seeder_save[3];
extern unsigned short int rng2_seeder[3];
extern unsigned short int rng2_seeder_save[3];

inline double gamma_ln(const double xx) {

#if __USE_PRECISE_GAMMA__==1
	static const double cof[14] = { 57.1562356658629235,
		-59.5979603554754912,
		14.1360979747417471,
		-0.491913816097620199,
		0.339946499848118887e-4,
		0.465236289270485756e-4,
		-0.983744753048795646e-4,
		0.158088703224912494e-3,
		-0.210264441724104883e-3,
		0.217439618115212643e-3,
		-0.164318106536763890e-3,
		0.844182239838527433e-4,
		-0.261908384015814087e-4,
		0.368991826595316234e-5
	};
#else
	static const double cof[6] = { 76.18009172947146,
		-86.50532032941677,
		24.01409824083091,
		-1.231739572450155,
		0.1208650973866179e-2,
		-0.5395239384953e-5 };
#endif

	double x, tmp, y, ser;
	int j;

	y = x = xx;

#if __USE_PRECISE_GAMMA__ == 1
	tmp = x + 5.24218750000000000;
	tmp = (x + 0.5) * log(tmp) - tmp;
	ser = 0.999999999999997092;
	for (j = 0;j < 14;j++) ser += cof[j] / ++y;
	return tmp + log(2.5066282746310005 * ser / x);

#else

	tmp = x + 5.5;
	tmp -= (x + 0.5) * log(tmp);
	ser = 1.000000000190015;
	for (j = 0; j <= 5; j++) ser += cof[j] / ++y;
	return -tmp + log(2.5066282746310005 * ser / x);

#endif

}

inline double sample_NormalSampler_0_1_0(void){
	double u, v, x, y, q;
	do {
		u = sample_uniform_rng2();
		v = 1.7156 * (sample_uniform_rng2() - 0.5);
		x = u - 0.449871;
		y = abs(v) + 0.386595;
		q = (x * x) + y * (0.19600 * y - 0.25472 * x);
	} while ((q > 0.27597) && (q > 0.27846 || (v * v) > -4.0 * log(u) * (u * u)));
	return(v / u);
}

// @summary Generate normal deviates using Ratio-of-Uniforms method
// @note Original implementation based on the method described in Numerical Recipes (3rd ed) 7.3.9
struct NormalSampler {
	double mu;
	double sigma;
	double sample(void) {

		double u, v, x, y, q;
		do {
			u = sample_uniform_rng2();
			v = 1.7156 * (sample_uniform_rng2() - 0.5);
			x = u - 0.449871;
			y = abs(v) + 0.386595;
			q = (x * x) + y * (0.19600 * y - 0.25472 * x);
		} while ((q > 0.27597) && (q > 0.27846 || (v * v) > -4.0 * log(u) * (u * u)));

		return(this->mu + this->sigma * v / u);
	}

};


inline NormalSampler* NormalSampler_init(const double mean, const double var, const int seed) {

    NormalSampler* norm = (NormalSampler*)malloc(sizeof(NormalSampler));

    norm->mu = mean;
    norm->sigma = var;

    return(norm);
}

inline void NormalSampler_destroy(NormalSampler* norm) {
    free(norm);
    norm = NULL;
}


// @brief Gamma1Sampler is a special case of GammaSampler with rate=1
// @note  Gamma1Sampler is used in BetaSampler
struct Gamma1Sampler {
	double alpha; // shape param

	bool alpha_changed; // true if alpha<1.0; false otherwise
	double old_alpha;   // set if alpha_changed==true

	// NormalSampler* normalSampler;

	double a1;
	double a2;

	double sample(void){
		double u, v, x, xsq;
		do {
			do {
				// x = this->normalSampler->sample();
				x = sample_NormalSampler_0_1_0();
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
};

inline Gamma1Sampler* Gamma1Sampler_init(const double shape) {

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
    // gamma->normalSampler = NormalSampler_init(0.0, 1.0, 0);

    return(gamma);
}

inline void Gamma1Sampler_destroy(Gamma1Sampler* gamma) {
    // NormalSampler_destroy(gamma->normalSampler);
    free(gamma);
    gamma = NULL;
}


// @brief Generate gamma deviates using Marsaglia-Tsang method
// @note 
struct GammaSampler {
	double alpha; // shape param
	double beta;  // rate param

	bool alpha_changed; // true if alpha<1.0; false otherwise
	double old_alpha;   // set if alpha_changed==true

	// NormalSampler* normalSampler;

	double a1;
	double a2;

	double sample(void){
		double u, v, x, xsq;
		do {
			do {
				// x = this->normalSampler->sample();
				x = sample_NormalSampler_0_1_0();
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

};

inline GammaSampler* GammaSampler_init(const double shape, const double rate) {

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
    // gamma->normalSampler = NormalSampler_init(0.0, 1.0, 0);

    return(gamma);
}

inline void GammaSampler_destroy(GammaSampler* gamma) {
    // NormalSampler_destroy(gamma->normalSampler);
    free(gamma);
    gamma = NULL;
}



struct PoissonSampler {

	double lm; // lambda
	double sq;
	double alxm;
	double g;
	bool st12; // lm < 12.0

};

inline PoissonSampler* PoissonSampler_init(const double lambda) {

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



inline void poissonSampler_sample_depths_same_mean(PoissonSampler* pois, int* n_sim_reads_arr, const int nSamples) {
	double em;
	double t;
	if (pois->st12) {

		for (int i = 0;i < nSamples;++i) {
			em = -1.0;
			t = 1.0;
			do {
				++em;
				t *= erand48(rng1_seeder);
			} while (t > (pois->g));
			n_sim_reads_arr[i] = em;
		}

	} else {
		double y;
		for (int i = 0;i < nSamples;++i) {
			do {
				do {
					y = tan(PI * erand48(rng1_seeder));
					em = (pois->sq) * y + (pois->lm);
				} while (em < 0.0);
				em = floor(em);
				t = 0.9 * (1.0 + y * y) * exp(em * (pois->alxm) - gamma_ln(em + 1.0) - (pois->g));
			} while (erand48(rng1_seeder) > t);

			n_sim_reads_arr[i] = em;
		}

	}
	return;
}

inline void poissonSampler_sample_depths_perSample_means(PoissonSampler** multipois, int* n_sim_reads_arr, const int nSamples) {
	double em;
	double t;
	PoissonSampler* pois = NULL;

	for (int i = 0;i < nSamples;++i) {
		pois = multipois[i];

		if (pois->st12) {
			em = -1.0;
			t = 1.0;
			do {
				++em;
				t *= erand48(rng1_seeder);
			} while (t > (pois->g));
			n_sim_reads_arr[i] = em;

		} else {
			double y;
			do {
				do {
					y = tan(PI * erand48(rng1_seeder));
					em = (pois->sq) * y + (pois->lm);
				} while (em < 0.0);
				em = floor(em);
				t = 0.9 * (1.0 + y * y) * exp(em * (pois->alxm) - gamma_ln(em + 1.0) - (pois->g));
			} while (erand48(rng1_seeder) > t);

			n_sim_reads_arr[i] = em;

		}
	}
	return;
}

#if __USE_STD_BETA__==1
struct BetaSampler {
	double alpha, beta;

	std::mt19937 generator;

	std::gamma_distribution<double>* gamma_x;
	std::gamma_distribution<double>* gamma_y;

	BetaSampler(const double mean, const double var, const int seed, FILE* arg_fp){

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

		fprintf(arg_fp, "\n-> Beta distribution shape parameters are estimated as alpha=%f and beta=%f (mean=%f, variance=%e)\n", this->alpha, this->beta, mean, var);

		this->gamma_x = new std::gamma_distribution<double>(this->alpha, 1.0);
		this->gamma_y = new std::gamma_distribution<double>(this->beta, 1.0);

		generator.seed(seed);
	}

	~BetaSampler() {
		delete gamma_x;
		delete gamma_y;
	}

	double sample(void) {
		std::gamma_distribution<double> gamma_x(alpha, 1.0);
		std::gamma_distribution<double> gamma_y(beta, 1.0);
		double x = gamma_x(generator);
		double y = gamma_y(generator);

		double ret = (x / (x + y));
		ASSERT(ret >= 0.0);
		ASSERT(ret <= 1.0);

		return (ret);
	}

};
#else
// @brief Generate beta deviates 
// @note Original implementation based on the method described in Numerical Recipes (3rd ed) 7.3.33
// @note Uses Gamma1Sampler
struct BetaSampler {
	double alpha;
	double beta;

	Gamma1Sampler* gamma_x;
	Gamma1Sampler* gamma_y;

	double sample(void){

		double x = this->gamma_x->sample();
		double y = this->gamma_y->sample();

		double ret = (x / (x + y));

		DEVASSERT(ret >= 0.0);
		DEVASSERT(ret <= 1.0);

		return (ret);
	}

};


inline BetaSampler* BetaSampler_init(const double mean, const double var, const int seed, FILE* arg_fp) {

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

    fprintf(arg_fp, "\n-> Beta distribution shape parameters are estimated as alpha=%f and beta=%f (mean=%f, variance=%e)\n", beta->alpha, beta->beta, mean, var);


    // beta->gamma_x = GammaSampler_init(beta->alpha, 1.0);
    // beta->gamma_y = GammaSampler_init(beta->beta, 1.0);
    beta->gamma_x = Gamma1Sampler_init(beta->alpha);
    beta->gamma_y = Gamma1Sampler_init(beta->beta);

    return(beta);
}

inline void BetaSampler_destroy(BetaSampler* beta) {
    Gamma1Sampler_destroy(beta->gamma_x);
    Gamma1Sampler_destroy(beta->gamma_y);
    free(beta);
    beta = NULL;
}


#endif


#endif // __RANDOM_GENERATOR__
