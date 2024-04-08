#ifndef __RANDOM_GENERATOR__
#define __RANDOM_GENERATOR__

#include <random> // std::mt19937, std::gamma_distribution


#include "shared.h"

#define sample_uniform_rng0() drand48()
#define sample_uniform_rng1() erand48(rng1_seeder)
#define sample_uniform_rng2() erand48(rng2_seeder)
#define sample_uniform_rng_rand ((double)rand() / (double)RAND_MAX + 1.0)
#define sample_from_range_rng_rand(min,max) ( (min) + (rand()) / (RAND_MAX / ((max) - (min) + 1) +1));

#define SQUARE(x) (((x)*(x)))


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

NormalSampler* NormalSampler_init(const double mu, const double sigma);
void NormalSampler_destroy(NormalSampler* norm);

// @brief Gamma1Sampler is a special case of GammaSampler with rate=1
// @note  Gamma1Sampler is used in BetaSampler
struct Gamma1Sampler {
	double alpha; // shape param

	bool alpha_changed; // true if alpha<1.0; false otherwise
	double old_alpha;   // set if alpha_changed==true

	NormalSampler* normalSampler;

	double a1;
	double a2;

	double sample();
};

Gamma1Sampler* Gamma1Sampler_init(const double shape);
void Gamma1Sampler_destroy(Gamma1Sampler* gamma);


// @brief Generate gamma deviates using Marsaglia-Tsang method
// @note 
struct GammaSampler {
	double alpha; // shape param
	double beta;  // rate param

	bool alpha_changed; // true if alpha<1.0; false otherwise
	double old_alpha;   // set if alpha_changed==true

	NormalSampler* normalSampler;

	double a1;
	double a2;

	double sample();
};

GammaSampler* GammaSampler_init(const double shape, const double rate);
void GammaSampler_destroy(GammaSampler* gamma);

struct PoissonSampler {

	double lm; // lambda
	double sq;
	double alxm;
	double g;
	bool st12; // lm < 12.0

	double sample(void);
	double sample_st12(void);
	double sample_bteq12(void);

};

PoissonSampler* PoissonSampler_init(const double lambda);

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

	BetaSampler(const double mean, const double var, const int seed);
	~BetaSampler();

	double sample();
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

	double sample();
};

//TODO rm seed?
BetaSampler* BetaSampler_init(const double mean, const double var, const int seed);
void BetaSampler_destroy(BetaSampler* beta);
#endif



#endif // __RANDOM_GENERATOR__
