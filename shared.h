#ifndef __SHARED__
#define __SHARED__

#include <stdio.h>   // fprintf
#include <stdlib.h>  // exit
#include <limits>    // std::numeric_limits

#include <float.h>  // DBL_MANT_DIG

#include "dev.h"

/* -> CONSTANTS --------------------------------------------------------------*/

const double NEG_INF = -std::numeric_limits<double>::infinity();

#define ARGS_DEPTH_ISUNDEF -1.0
#define ARGS_DEPTH_ISFILE -2.0
#define ARGS_DEPTH_ISINF -3.0

// args->useUnknownAllele
// notations for representing non-reference unobserved alleles
// 	ARGS_NONREF_NOTATION_STAR		<*> symbolic alternate allele (used in bcftools)
// 	ARGS_NONREF_NOTATION_NON_REF	<NON_REF> gVCF NON_REF notation (used in GATK)
#define ARGS_NONREF_NOTATION_EXPLODE_ACGT 0
#define ARGS_NONREF_NOTATION_STAR 1
#define ARGS_NONREF_NOTATION_NON_REF 2


// precalculated value for log10(3)
#define PRE_CALC_LOG10_3 0.47712125471966244

#define MAXGL 0
#define MINGL NEG_INF

#define MAXPL 255
#define MINPL 0

#define MAXGP 1.0
#define MINGP 0.0

#define SIM_FORWARD_STRAND 0
#define SIM_REVERSE_STRAND 1

// maximum number of genotypes to simulate
#define MAX_NGTS 10  // {AA,AC,CC,AG,CG,GG,AT,CT,GT,TT}

// maximum number of alleles to simulate
#define MAX_NALLELES 4  // {A,C,G,T}

// assume diploid samples
#define SIM_PLOIDY 2

// source: bcftools/bam2bcf.c L41
#define CAP_DIST 25

// source: bcftools/bam2bcf.c L381
#define CAP_BASEQ 63

#define BCF_GT_PHASED_0 3  // bcf_gt_phased(0)
#define BCF_GT_PHASED_1 5  // bcf_gt_phased(1)

/* -> FUNCTION-LIKE MACROS ---------------------------------------------------*/

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

/*
 * Macro:[DBL_MAX_DIG_TOPRINT]
 * 	maximum number of digits needed to print a double
 *
 * 	longest number == smalles negative number
 * 		-pow(2, DBL_MIN_EXP - DBL_MANT_DIG)
 * 	-pow(2,-N) needs 3+N digits
 * 		to represent (sign, decimal point, N digits)
 * 		'-0.<N digits>'
 *
 * @requires <float.h>
 */
#define DBL_MAX_DIG_TOPRINT 3 + DBL_MANT_DIG - DBL_MIN_EXP

#define DBL_MAXDIG10 (2 + (DBL_MANT_DIG * 30103UL) / 100000UL)
/*
 * Macro:[AT]
 * inject the file and line info as string
 */
#define STRINGIFY(x) #x
#define ASSTR(x) STRINGIFY(x)
#define AT __FILE__ ":" ASSTR(__LINE__)

/*
 * Macro:[ERROR]
 * print a custom error message and exit the program
 */
#define ERROR(...)                                                           \
    do {                                                                     \
        fprintf(stderr, "\n\n*******\n[ERROR](%s)<%s:%d>\n\t", __FUNCTION__, \
                __FILE__, __LINE__);                                         \
        fprintf(stderr, __VA_ARGS__);                                        \
        fprintf(stderr, "\n*******\n");                                      \
        exit(1);                                                             \
    } while (0);

/*
 * Macro:[NEVER]
 * indicates that a point in the code should never be reached
 */
#define NEVER                                                               \
    do {                                                                    \
        ERROR(                                                              \
            "Control should never reach this point; please report this to " \
            "the developers.")                                              \
    } while (0);

/*
 * Macro:[ASSERT]
 * evaluate an expression, works the same way as the C-macro assert
 * except that DEBUG does not affect it (it is always active)
 * also prints the file and line info and exits the program
 * if the expression evaluates to false
 */
#define ASSERT(expr)                                                        \
    do {                                                                    \
        if (!((expr))) {                                                    \
            fprintf(stderr, "\n\n*******\n[ERROR](%s/%s:%d) %s\n*******\n", \
                    __FILE__, __FUNCTION__, __LINE__, #expr);               \
            exit(1);                                                        \
        }                                                                   \
    } while (0);

/*
 * Macro:[WARN]
 * print a custom warning message
 */
#define WARN(...)                                                            \
    do {                                                                     \
        fprintf(stderr, "\n\n[WARNING](%s/%s:%d): ", __FILE__, __FUNCTION__, \
                __LINE__);                                                   \
        fprintf(stderr, __VA_ARGS__);                                        \
        fprintf(stderr, "\n");                                               \
    } while (0);

/*
 * Macro:[VWARN]
 * print a custom warning message if verbose is set
 */
#define VWARN(...)                                                 \
    do {                                                           \
        if (0 != args->verbose) {                                  \
            fprintf(stderr, "\n\n[WARNING](%s/%s:%d): ", __FILE__, \
                    __FUNCTION__, __LINE__);                       \
            fprintf(stderr, __VA_ARGS__);                          \
            fprintf(stderr, "\n");                                 \
        }                                                          \
    } while (0);

/* LOOKUP TABLES & LOOKUP FUNCTIONS ------------------------------------------*/

// @brief lut_qs_to_qs2 - map quality score to squared quality score
extern const int lut_qs_to_qs2[64];

// [R]
// qToP<-function(q){10^(-q/10)}
// CAP_BASEQ=63
// arrSize=CAP_BASEQ+1
// paste0("qScore_to_errorProb[",arrSize,"]={",paste(unlist(format(lapply(0:CAP_BASEQ,FUN=qToP),scientific=F)),collapse=","),"};")
// [1]
// "qScore_to_errorProb[64]={1,0.7943282,0.6309573,0.5011872,0.3981072,0.3162278,0.2511886,0.1995262,0.1584893,0.1258925,0.1,0.07943282,0.06309573,0.05011872,0.03981072,0.03162278,0.02511886,0.01995262,0.01584893,0.01258925,0.01,0.007943282,0.006309573,0.005011872,0.003981072,0.003162278,0.002511886,0.001995262,0.001584893,0.001258925,0.001,0.0007943282,0.0006309573,0.0005011872,0.0003981072,0.0003162278,0.0002511886,0.0001995262,0.0001584893,0.0001258925,0.0001,0.00007943282,0.00006309573,0.00005011872,0.00003981072,0.00003162278,0.00002511886,0.00001995262,0.00001584893,0.00001258925,0.00001,0.000007943282,0.000006309573,0.000005011872,0.000003981072,0.000003162278,0.000002511886,0.000001995262,0.000001584893,0.000001258925,0.000001,0.0000007943282,0.0000006309573,0.0000005011872};"
extern const double qScore_to_errorProb[64];

// qScore	phred-scaled quality score
// 			qScore = -10 * log10(error_probability)

// @brief qs_to_qs2 - get squared quality score from quality score
// @details
// [R]
// > n <- 63
// > paste0("const int qs_to_qs2[", n + 1, "] = {", paste(sapply(0:n,
// function(n) n*n),collapse = ", "),"};") [1] "const int qs_to_qs2[64] = {0, 1,
// 4, 9, 16, 25, 36, 49, 64, 81, 100, 121, 144, 169, 196, 225, 256, 289, 324,
// 361, 400, 441, 484, 529, 576, 625, 676, 729, 784, 841, 900, 961, 1024, 1089,
// 1156, 1225, 1296, 1369, 1444, 1521, 1600, 1681, 1764, 1849, 1936, 2025, 2116,
// 2209, 2304, 2401, 2500, 2601, 2704, 2809, 2916, 3025, 3136, 3249, 3364, 3481,
// 3600, 3721, 3844, 3969};"
inline int qs_to_qs2(const int q) {
    if (0 == q) {
        return 0;
    } else if (q > CAP_BASEQ) {
        return lut_qs_to_qs2[CAP_BASEQ];
    } else if (q > 0) {
        return lut_qs_to_qs2[q];
    } else {
        NEVER;
    }
}

// @brief lut_myGtIdx_to_vcfGtIdx
//      map internal representation order to vcf genotype order
// 		for 4 alleles {A,C,G,T}
//
// FROM: Internal representation
//
// | 0 | 1 | 2 | 3 |
// | A | C | G | T |
//
// | 0  | 1  | 2  | 3  | 4  | 5  | 6  | 7  | 8  | 9  |
// | AA | AC | AG | AT | CC | CG | GG | CT | GT | TT |
//
// TO: VCF genotype order
//
// for P=ploidy and N=number of alternate alleles;
// [pseudocode] for a_p in 0:N; for a_p-1 in 0:a_p; print(a1,a2);
//
// For P=2 N=3
//
// | 0  | 1  | 4  | 2  | 5  | 7  | 3  | 6  | 8   | 9   |
// | AA | AC | CC | AG | CG | GG | AT | CT | GT  | TT  |
//
extern const int lut_myGtIdx_to_vcfGtIdx[10];

// @brief myGtIdx_to_vcfGtIdx - get vcf genotype index from internal
// representation genotype index
// @param 	gti (int) - genotype index in internal representation order
// (0-9:AA,AC,AG,AT,CC,CG,GG,CT,GT,TT)
// @return	    (int) - genotype index in vcf order
// (0-9:AA,AC,CC,AG,CG,GG,AT,CT,GT,TT)
inline int myGtIdx_to_vcfGtIdx(const int gti) {
    ASSERT(gti >= 0 && gti < 10);
    return (lut_myGtIdx_to_vcfGtIdx[gti]);
}

// @brief myAllele_to_vcfGtIdx - get vcf genotype index from internal
// representation alleles
// @param 	a1 (int) - allele 1 in internal representation order
// (0-3:A,C,G,T)
// @param 	a2 (int) - allele 2 in internal representation order
// (0-3:A,C,G,T)
// @return	   (int) - genotype index in vcf order
// (0-9:AA,AC,CC,AG,CG,GG,AT,CT,GT,TT)
inline int myAllele_to_vcfGtIdx(const int a1, const int a2) {
    if (a1 > a2)
        return (lut_myGtIdx_to_vcfGtIdx[a2 * 4 + a1]);
    // else ---> if (a1 <= a2)
    return (lut_myGtIdx_to_vcfGtIdx[a1 * 4 + a2]);
}

// @brief nAlleles_to_nGenotypes - map number of alleles to number of genotypes
// @details
//      for (i=0;i<5;++i)
//          nAlleles_to_nGenotypes[i] == nAlleles2nGenotypes(i) == i*(i+1)/2
extern const int lut_nAlleles_to_nGenotypes[6];

// @brief nAlleles_to_nGenotypes - get number of genotypes from number of
// alleles
// @param 	n (int) - number of alleles
// @return	  (int) - number of unique unordered genotypes expected for n
// alleles
//                    e.g. A, C -> AA AC CC -> 3
inline int nAlleles_to_nGenotypes(const int n) {
    ASSERT(n > 0 && n < 6);
    return (lut_nAlleles_to_nGenotypes[n]);
}

// @brief lut_acgt_offsets - map alleles to genotype offsets (internal
// representation)
// @details
//
// | 0  | 1  | 2  | 3  | 4  | 5  | 6  | 7  | 8  | 9  |
// | AA | AC | AG | AT | CC | CG | GG | CT | GT | TT |
//
// Matrix format
// |    | 0  | 1  | 2  | 3  |		| 4  | 5  | 6  | 7  | 8  | 9  |
// |----|----|----|----|----|		|----|----|----|----|----|----|
// | 0  | AA | AC | AG | AT |		| CC | CG | GG | CT | GT | TT |
// | 1  | CC | AC | CG | CT |		| AA | AG | AT | GG | GT | TT |
// | 2  | GG | AG | CG | GT |		| AA | AC | AT | CC | CT | TT |
// | 3  | TT | AT | CT | GT | 		| AA | AC | AG | CC | CG | GG |
//
extern const int lut_acgt_offsets[4][10];

#endif  // __SHARED__
