#ifndef __GL_METHODS__
#define __GL_METHODS__

#include "bcf_utils.h"

typedef struct glModel1Struct glModel1Struct;

struct glModel1Struct {
    errmod_t* errmod;
};
extern glModel1Struct* glModel1;
glModel1Struct* glModel1_init(void);
void glModel1_destroy(glModel1Struct* glModel1);

void alleles_calculate_gls_log10_glModel2_fixedQScore(simRecord* sim);
void alleles_calculate_gls_log10_glModel2_precise0(simRecord* sim);
void alleles_calculate_gls_log10_glModel2_precise1(simRecord* sim);
void alleles_calculate_gls_log10_glModel1_fixedQScore(simRecord* sim);
void alleles_calculate_gls_log10_glModel1(simRecord* sim);

#endif  // __GL_METHODS__
