#include "gl_methods.h"
#include "io.h"

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

    int** qScores = NULL;
    qScores = (PROGRAM_WILL_ADJUST_QS_FOR_GL) ? sim->adj_base_qScores : sim->base_qScores;

    for (int s = 0;s < sim->nSamples;s++) {
        n_sim_reads = sim->fmt_dp_arr[s];

        sample_gl_arr = sim->gl_arr + (s * sim->nGenotypes);

        if (0 != n_sim_reads) {

            for (int read_i = 0; read_i < n_sim_reads; read_i++) {

                bo = sim->bases[s][read_i];

                qs = qScores[s][read_i];

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

    int** qScores = NULL;
    qScores = (PROGRAM_WILL_ADJUST_QS_FOR_GL) ? sim->adj_base_qScores : sim->base_qScores;

    // 5 * 5 = 25 full allele combinations
    // for 5 alleles (A,C,G,T,<*>) 
    float fpls[25];

    for (int s = 0;s < sim->nSamples;s++) {
        n_sim_reads = sim->fmt_dp_arr[s];

        sample_gl_arr = sim->gl_arr + (s * sim->nGenotypes);

        if (0 != n_sim_reads) {

            uint16_t ubases[n_sim_reads];

            for (int i = 0;i < n_sim_reads;++i) {
                qs = qScores[s][i];
                ubases[i] = qs << 5 | sim->bases[s][i];
            }

            // 5 alleles (A,C,G,T,<*>)
            errmod_cal(args->gl1errmod, n_sim_reads, 5, ubases, fpls);

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

    int q5 = (PROGRAM_WILL_ADJUST_QS_FOR_GL) ? args->preCalc->adj_q5 : args->preCalc->q5;

    for (int s = 0;s < sim->nSamples;s++) {
        n_sim_reads = sim->fmt_dp_arr[s];

        sample_gl_arr = sim->gl_arr + (s * sim->nGenotypes);

        if (0 != n_sim_reads) {

            uint16_t ubases[n_sim_reads];
            for (int i = 0;i < n_sim_reads;++i) {
                ubases[i] = q5 | sim->bases[s][i];
            }

            // 5 alleles (A,C,G,T,<*>)
            errmod_cal(args->gl1errmod, n_sim_reads, 5, ubases, fpls);

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
