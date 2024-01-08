#include "gl_methods.h"
#include "io.h"

extern glModel1Struct* glModel1;
glModel1Struct* glModel1_init(void) {
    glModel1Struct* glModel1 = (glModel1Struct*)malloc(sizeof(glModel1Struct));
    ASSERT(NULL != glModel1);

    glModel1->errmod = errmod_init(1.0 - args->glModel1_theta);
    ASSERT(NULL != glModel1->errmod);

    return(glModel1);
}

void glModel1_destroy(glModel1Struct* glModel1) {
    errmod_destroy(glModel1->errmod);
    free(glModel1);
    return;
}


void alleles_calculate_gls_log10_glModel2_fixedQScore(simRecord* sim) {

    const double homT = args->preCalc->homT;
    const double het = args->preCalc->het;
    const double homF = args->preCalc->homF;

    float* sample_gl_arr = NULL;
    int n_sim_reads = 0;
    int bo;
    int ao;

    for (int s = 0;s < sim->nSamples;s++) {
        n_sim_reads = sim->fmt_dp_arr[s];

        sample_gl_arr = sim->gl_arr + (s * sim->nGenotypes);

        if (0 != n_sim_reads) {

            for (int read_i = 0; read_i < n_sim_reads; read_i++) {

                bo = sim->bases[s][read_i];
                ao = sim->acgt2alleles[bo];

                sample_gl_arr[bcf_alleles2gt(ao, ao)] += homT;  // homozygotic hit
                for (int a1 = 0;a1 < sim->nAlleles;++a1) {
                    if (a1 == ao) {
                        continue;
                    } else if (a1 > ao) {
                        sample_gl_arr[bcf_alleles2gt(a1, ao)] += het;  // heterozygotic hit
                        for (int a2 = a1;a2 < sim->nAlleles;++a2) {
                            if (a2 != ao) {
                                sample_gl_arr[bcf_alleles2gt(a2, a1)] += homF;  // non hit
                            }
                        }

                    } else if (a1 < ao) {
                        sample_gl_arr[bcf_alleles2gt(ao, a1)] += het;  // heterozygotic hit

                        for (int a2 = a1;a2 < sim->nAlleles;++a2) {
                            if (a2 != ao) {
                                sample_gl_arr[bcf_alleles2gt(a2, a1)] += homF;  // non hit
                            }
                        }
                    }
                }

                float max = NEG_INF;
                for (int i = 0;i < sim->nGenotypes;++i) {
                    if (sample_gl_arr[i] > max) {
                        max = sample_gl_arr[i];
                    }
                }
                for (int i = 0;i < sim->nGenotypes;++i) {
                    sample_gl_arr[i] -= max;
                }
            }
        } else {
            // set all to missing for ind

            for (int g = 0;g < sim->nGenotypes;++g) {
                bcf_float_set_missing(sample_gl_arr[g]);
            }
        }
    }

}

void alleles_calculate_gls_log10_glModel2_precise0(simRecord* sim) {

    double homT;
    double het;
    double homF;

    float* sample_gl_arr = NULL;
    int n_sim_reads = 0;
    int bo;
    int ao;
    int qs;

    for (int s = 0;s < sim->nSamples;s++) {
        n_sim_reads = sim->fmt_dp_arr[s];

        sample_gl_arr = sim->gl_arr + (s * sim->nGenotypes);

        if (0 != n_sim_reads) {

            for (int read_i = 0; read_i < n_sim_reads; read_i++) {

                bo = sim->bases[s][read_i];

                qs = sim->base_qScores[s][read_i];
                homT = qScore_to_log10_gl[0][qs];
                het = qScore_to_log10_gl[1][qs];
                homF = qScore_to_log10_gl[2][qs];


                ao = sim->acgt2alleles[bo];
                sample_gl_arr[bcf_alleles2gt(ao, ao)] += homT;  // homozygotic hit
                for (int a1 = 0;a1 < sim->nAlleles;++a1) {
                    if (a1 == ao) {
                        continue;
                    } else if (a1 > ao) {
                        sample_gl_arr[bcf_alleles2gt(a1, ao)] += het;  // heterozygotic hit
                        for (int a2 = a1;a2 < sim->nAlleles;++a2) {
                            if (a2 != ao) {
                                sample_gl_arr[bcf_alleles2gt(a2, a1)] += homF;  // non hit
                            }
                        }

                    } else if (a1 < ao) {
                        sample_gl_arr[bcf_alleles2gt(ao, a1)] += het;  // heterozygotic hit

                        for (int a2 = a1;a2 < sim->nAlleles;++a2) {
                            if (a2 != ao) {
                                sample_gl_arr[bcf_alleles2gt(a2, a1)] += homF;  // non hit
                            }
                        }
                    }
                }

                float max = NEG_INF;
                for (int i = 0;i < sim->nGenotypes;++i) {
                    if (sample_gl_arr[i] > max) {
                        max = sample_gl_arr[i];
                    }
                }
                for (int i = 0;i < sim->nGenotypes;++i) {
                    sample_gl_arr[i] -= max;
                }
            }
        } else {
            // set all to missing for ind

            for (int g = 0;g < sim->nGenotypes;++g) {
                bcf_float_set_missing(sample_gl_arr[g]);
            }

        }
    }

}

void alleles_calculate_gls_log10_glModel2_precise1(simRecord* sim) {

    double homT;
    double het;
    double homF;
    double e;

    float* sample_gl_arr = NULL;
    int n_sim_reads = 0;
    int bo;
    int ao;

    for (int s = 0;s < sim->nSamples;s++) {
        n_sim_reads = sim->fmt_dp_arr[s];

        sample_gl_arr = sim->gl_arr + (s * sim->nGenotypes);

        if (0 != n_sim_reads) {

            for (int read_i = 0; read_i < n_sim_reads; read_i++) {

                bo = sim->bases[s][read_i];

                e = sim->base_error_probs[s][read_i];

                if (0.0 == e) {
                    homT = 0.0;
                    het = -0.30103;
                    homF = NEG_INF;
                } else {
                    homT = log10(1.0 - e);
                    het = log10((1.0 - e) / 2.0 + e / 6.0);
                    homF = log10(e / 3.0);
                }

                ao = sim->acgt2alleles[bo];
                sample_gl_arr[bcf_alleles2gt(ao, ao)] += homT;  // homozygotic hit
                for (int a1 = 0;a1 < sim->nAlleles;++a1) {
                    if (a1 == ao) {
                        continue;
                    } else if (a1 > ao) {
                        sample_gl_arr[bcf_alleles2gt(a1, ao)] += het;  // heterozygotic hit
                        for (int a2 = a1;a2 < sim->nAlleles;++a2) {
                            if (a2 != ao) {
                                sample_gl_arr[bcf_alleles2gt(a2, a1)] += homF;  // non hit
                            }
                        }

                    } else if (a1 < ao) {
                        sample_gl_arr[bcf_alleles2gt(ao, a1)] += het;  // heterozygotic hit

                        for (int a2 = a1;a2 < sim->nAlleles;++a2) {
                            if (a2 != ao) {
                                sample_gl_arr[bcf_alleles2gt(a2, a1)] += homF;  // non hit
                            }
                        }
                    }
                }

                float max = NEG_INF;
                for (int i = 0;i < sim->nGenotypes;++i) {
                    if (sample_gl_arr[i] > max) {
                        max = sample_gl_arr[i];
                    }
                }
                for (int i = 0;i < sim->nGenotypes;++i) {
                    sample_gl_arr[i] -= max;
                }
            }
        } else {
            // set all to missing for ind

            for (int g = 0;g < sim->nGenotypes;++g) {
                bcf_float_set_missing(sample_gl_arr[g]);
            }

        }
    }

}

void alleles_calculate_gls_log10_glModel1(simRecord* sim) {

    float max = NEG_INF;

    float* sample_gl_arr = NULL;
    int n_sim_reads = 0;
    int qs;

    int a1, a2, b1, b2;
    int g = 0;

    // 5 * 5 = 25 full allele combinations
    // for 5 alleles (A,C,G,T,<*>) 
    float fpls[25];

    for (int s = 0;s < sim->nSamples;s++) {
        n_sim_reads = sim->fmt_dp_arr[s];

        sample_gl_arr = sim->gl_arr + (s * sim->nGenotypes);

        if (0 != n_sim_reads) {

            uint16_t ubases[n_sim_reads];

            for (int i = 0;i < n_sim_reads;++i) {
                qs = sim->base_qScores[s][i];
                ubases[i] = qs << 5 | sim->bases[s][i];
            }

            // 5 alleles (A,C,G,T,<*>)
            errmod_cal(glModel1->errmod, n_sim_reads, 5, ubases, fpls);

            max = NEG_INF;
            g = 0;

            for (a2 = 0; a2 < sim->nAlleles; ++a2) {
                b2 = sim->alleles2acgt[a2];
                for (a1 = 0; a1 <= a2; ++a1) {
                    b1 = sim->alleles2acgt[a1];

                    sample_gl_arr[g] = ((-1.0 * fpls[(b1 * 5 + b2)]) / 10.0);

                    if (sample_gl_arr[g] > max) {
                        max = sample_gl_arr[g];
                    }

                    ++g;
                }
            }

            DEVASSERT(g == sim->nGenotypes);

            for (int i = 0;i < g;++i) {
                sample_gl_arr[i] -= max;
            }

        } else {
            // set all to missing for ind

            for (int g = 0;g < sim->nGenotypes;++g) {
                bcf_float_set_missing(sample_gl_arr[g]);
            }

        }
    }

}

void alleles_calculate_gls_log10_glModel1_fixedQScore(simRecord* sim) {

    float max = NEG_INF;

    float* sample_gl_arr = NULL;
    int n_sim_reads = 0;

    int a1, a2, b1, b2;
    int g = 0;

    // 5 * 5 = 25 full allele combinations
    // for 5 alleles (A,C,G,T,<*>) 
    float fpls[25];

    for (int s = 0;s < sim->nSamples;s++) {
        n_sim_reads = sim->fmt_dp_arr[s];

        sample_gl_arr = sim->gl_arr + (s * sim->nGenotypes);

        if (0 != n_sim_reads) {

            uint16_t ubases[n_sim_reads];
            for (int i = 0;i < n_sim_reads;++i) {
                ubases[i] = args->preCalc->q5 | sim->bases[s][i];
            }

            // 5 alleles (A,C,G,T,<*>)
            errmod_cal(glModel1->errmod, n_sim_reads, 5, ubases, fpls);

            max = NEG_INF;
            g = 0;

            for (a2 = 0; a2 < sim->nAlleles; ++a2) {
                b2 = sim->alleles2acgt[a2];
                for (a1 = 0; a1 <= a2; ++a1) {
                    b1 = sim->alleles2acgt[a1];

                    sample_gl_arr[g] = ((-1.0 * fpls[(b1 * 5 + b2)]) / 10.0);

                    if (sample_gl_arr[g] > max) {
                        max = sample_gl_arr[g];
                    }

                    ++g;
                }
            }

            DEVASSERT(g == sim->nGenotypes);

            for (int i = 0;i < g;++i) {
                sample_gl_arr[i] -= max;
            }

        } else {
            // set all to missing for ind

            for (int g = 0;g < sim->nGenotypes;++g) {
                bcf_float_set_missing(sample_gl_arr[g]);
            }

        }
    }

}


// --- BELOW : WITH OLD GL_VALS ARRAY --

// glModel2 using non-precise gl error probability using qScore to gl lut
void calculate_gls_log10_glModel2_precise0(const int bo, const int qScore, double* sample_acgt_gls) {

    double homT = qScore_to_log10_gl[0][qScore];
    double het = qScore_to_log10_gl[1][qScore];
    double homF = qScore_to_log10_gl[2][qScore];

#if DEV==1
    int n = 0; //DEBUG
#endif

    // bcf_alelles2gt(bo, bo) = ((bo * (bo + 1) / (2 + bo));
    sample_acgt_gls[(bo * (bo + 1) / 2 + bo)] += homT;  // homozygotic hit

#if DEV==1
    n++;
#endif


    for (int b1 = 0;b1 < 5;++b1) {
        if (b1 == bo) {
            continue;
        } else if (b1 > bo) {
            // if b1>bo; then bcf_alleles2gt(b1,bo) == b1 * (b1 + 1) / 2 + bo
            sample_acgt_gls[(b1 * (b1 + 1) / 2 + bo)] += het;  // heterozygotic hit

#if DEV==1
            n++;
#endif

            for (int b2 = b1;b2 < 5;++b2) {
                // always: b2>=b1
                if (b2 != bo) {
                    // if b2>b1; then bcf_alleles2gt(b1,b2) == b2 * (b2 + 1) / 2 + b1
                    sample_acgt_gls[(b2 * (b2 + 1) / 2 + b1)] += homF;  // non hit

#if DEV==1
                    n++;
#endif

                }
            }

        } else if (b1 < bo) {
            // if b1<bo; then bcf_alleles2gt(b1,bo) == bo * (bo + 1) / 2 + b1
            sample_acgt_gls[(bo * (bo + 1) / 2 + b1)] += het;  // heterozygotic hit

#if DEV==1
            n++;
#endif

            for (int b2 = b1;b2 < 5;++b2) {
                // always: b2>=b1
                if (b2 != bo) {
                    // if b2>b1; then bcf_alleles2gt(b1,b2) == b2 * (b2 + 1) / 2 + b1
                    sample_acgt_gls[(b2 * (b2 + 1) / 2 + b1)] += homF;  // non hit

#if DEV==1
                    n++;
#endif

                }
            }
        }
    }

#if DEV==1
    ASSERT(n == 15);
#endif

}


// glModel2 using precise gl error probability from const double e
void calculate_gls_log10_glModel2_precise1(const int bo, const double e, double* sample_acgt_gls) {

    double homT = 0.0;
    double het = 0.0;
    double homF = 0.0;

    if (0.0 == e) {
        homT = 0.0;
        het = -0.30103;
        homF = NEG_INF;
    } else {
        homT = log10(1.0 - e);
        het = log10((1.0 - e) / 2.0 + e / 6.0);
        homF = log10(e / 3.0);
    }

#if DEV==1
    int n = 0; //DEBUG
#endif


    // bcf_alelles2gt(bo, bo) = ((bo * (bo + 1) / (2 + bo));
    sample_acgt_gls[(bo * (bo + 1) / 2 + bo)] += homT;  // homozygotic hit

#if DEV==1
    n++;
#endif


    for (int b1 = 0;b1 < 5;++b1) {
        if (b1 == bo) {
            continue;
        } else if (b1 > bo) {
            // if b1>bo; then bcf_alleles2gt(b1,bo) == b1 * (b1 + 1) / 2 + bo
            sample_acgt_gls[(b1 * (b1 + 1) / 2 + bo)] += het;  // heterozygotic hit

#if DEV==1
            n++;
#endif


            for (int b2 = b1;b2 < 5;++b2) {
                // always: b2>=b1
                if (b2 != bo) {
                    // if b2>b1; then bcf_alleles2gt(b1,b2) == b2 * (b2 + 1) / 2 + b1
                    sample_acgt_gls[(b2 * (b2 + 1) / 2 + b1)] += homF;  // non hit

#if DEV==1
                    n++;
#endif

                }
            }

        } else if (b1 < bo) {
            // if b1<bo; then bcf_alleles2gt(b1,bo) == bo * (bo + 1) / 2 + b1
            sample_acgt_gls[(bo * (bo + 1) / 2 + b1)] += het;  // heterozygotic hit

#if DEV==1
            n++;
#endif


            for (int b2 = b1;b2 < 5;++b2) {
                // always: b2>=b1
                if (b2 != bo) {
                    // if b2>b1; then bcf_alleles2gt(b1,b2) == b2 * (b2 + 1) / 2 + b1
                    sample_acgt_gls[(b2 * (b2 + 1) / 2 + b1)] += homF;  // non hit

#if DEV==1
                    n++;
#endif

                }
            }
        }
    }

#if DEV==1
    ASSERT(n == 15);
#endif

}

// glModel1 with fixed qScore
void calculate_gls_log10_glModel1(int* bases, int n_bases, float* gls, int* qScores, glModel1Struct* glModel1, simRecord* sim) {

    uint16_t ubases[n_bases];

    // 5 * 5 = 25 full allele combinations
    // for 5 alleles (A,C,G,T,<*>) 
    float fpls[25];

    int q;

    for (int i = 0;i < n_bases;++i) {
        q = qScores[i];
        ubases[i] = q << 5 | bases[i];
    }

    // 5 alleles (A,C,G,T,<*>)
    errmod_cal(glModel1->errmod, n_bases, 5, ubases, fpls);

    float max = NEG_INF;
    int a1, a2;
    int b1, b2;
    int g = 0;

    for (a2 = 0; a2 < sim->nAlleles; ++a2) {
        b2 = sim->alleles2acgt[a2];
        for (a1 = 0; a1 <= a2; ++a1) {
            b1 = sim->alleles2acgt[a1];

            gls[g] = ((-1.0 * fpls[(b1 * 5 + b2)]) / 10.0);

            if (gls[g] > max) {
                max = gls[g];
            }

            ++g;
        }
    }

    DEVASSERT(g == sim->nGenotypes);

    for (int i = 0;i < g;++i) {
        gls[i] -= max;
    }

}

