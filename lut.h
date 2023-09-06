#ifndef __LUT__
#define __LUT__

// lookup table for squared quality scores 
extern const int lut_qscore2[64];


/*
 * Macro:[ QSCORE2(q) ]
 * Get the squared quality score
 */
#define QSCORE2(q) ( ( (q) > 64 ) ? 4096 : lut_qscore2[(q)] )


extern const char SIM_ALLELES[5];

extern const int acgt_genotype_order_lut[10];

extern const int offsets[4][10];


#endif // __LUT__
