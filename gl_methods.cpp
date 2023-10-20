#include "gl_methods.h"
#include "io.h"

// method1: direct genotype likelihood method
void calc_gl_log10_method1_precalc(int base, double e, double* acgt_gls) {
    /*0=AA, 1=AC, 2=AG, 3=AT, 4=CC, 5=CG, 6=CT, 7=GG, 8= GT, 9= TT*/

    /*0=AA*/
    acgt_gls[lut_acgt_offsets[base][0]] += args->homT;  // homozygotic hit

    /*1=AC, 2=AG, 3=AT*/
    for (int o = 1; o < 4; o++)  // heterozygotic hit
        acgt_gls[lut_acgt_offsets[base][o]] += args->het;

    /*4=CC, 5=CG, 6=CT, 7=GG, 8= GT, 9= TT*/
    for (int o = 4; o < 10; o++)  // non hit
        acgt_gls[lut_acgt_offsets[base][o]] += args->homF;
}

// method1: direct genotype likelihood method
void calc_gl_log10_method1(int base, double e, double* acgt_gls) {
    /*0=AA, 1=AC, 2=AG, 3=AT, 4=CC, 5=CG, 6=CT, 7=GG, 8= GT, 9= TT*/

    double homT = log10(1.0 - e);
    double het = log10((1.0 - e) / 2.0 + e / 6.0);
    double homF = log10(e / 3.0);

    /*0=AA*/
    acgt_gls[lut_acgt_offsets[base][0]] += homT;  // homozygotic hit

    /*1=AC, 2=AG, 3=AT*/
    for (int o = 1; o < 4; o++)  // heterozygotic hit
        acgt_gls[lut_acgt_offsets[base][o]] += het;

    /*4=CC, 5=CG, 6=CT, 7=GG, 8= GT, 9= TT*/
    for (int o = 4; o < 10; o++)  // non hit
        acgt_gls[lut_acgt_offsets[base][o]] += homF;
}

void rescale_likelihood_ratio(double* like) {
    // rescale to likeratios
    double mx = like[0];

    for (int i = 1; i < 10; i++) {
        if (like[i] > mx) {
            mx = like[i];
        }
    }

    for (int i = 0; i < 10; i++) {
        like[i] -= mx;
    }
}

void rescale_likelihood_ratio(float* like) {
    // rescale to likeratios
    float mx = like[0];

    for (int i = 1; i < 10; i++) {
        if (like[i] > mx) {
            mx = like[i];
        }
    }

    for (int i = 0; i < 10; i++) {
        like[i] -= mx;
    }
}

void rescale_likelihood_ratio(float* like, const int size) {
    ASSERT(size > 0 && size <= MAX_NGTS);
    float mx = like[0];

    for (int i = 1; i < size; i++) {
        if (like[i] > mx) {
            mx = like[i];
        }
    }

    for (int i = 0; i < size; i++) {
        like[i] -= mx;
    }
}
