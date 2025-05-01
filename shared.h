#ifndef __SHARED__
#define __SHARED__

#include <stdio.h>   /// fprintf
#include <stdlib.h>  /// exit
#include <limits>    /// std::numeric_limits
#include <float.h>  /// DBL_MANT_DIG

#include "dev.h"


/* -> CONSTANTS --------------------------------------------------------------*/

/// --> CONSTANTS:RNG

/// source: rand48/rand48.h
#define VCFGL_RAND48_SEED_0   (0x330e)
#define VCFGL_RAND48_SEED_1   (0xabcd)
#define VCFGL_RAND48_SEED_2   (0x1234)

/// init to { VCFGL_RAND48_SEED_0, VCFGL_RAND48_SEED_1, VCFGL_RAND48_SEED_2 };
#define SEEDER_INIT { 0x330e, 0xabcd, 0x1234 }

#define KS_INIT { 0, 0, NULL } 


#ifndef __USE_PRECISE_GAMMA__
#define __USE_PRECISE_GAMMA__ 0
#endif

#ifndef __USE_STD_BETA__
#define __USE_STD_BETA__ 1
#endif

/// --> CONSTANTS:MATH

#define PI 3.141592654

const double NEG_INF = -std::numeric_limits<double>::infinity();

#define PRE_CALC_LOG10_3 (0.47712125471966244)
#define PRE_CALC_LN_3 (1.0986122886681098)

#define QSCORE_PHRED_ENCODING_OFFSET 33

/// --> CONSTANTS:BUFFER SIZE

#define KSTRING_BGZF_WRITE_BUFFER_SIZE 4096

#define BUFSIZE_NBASES 50

#define BUFSIZE_NINDS 100


/// --> CONSTANTS:ARGS

#define ARG_DEPTH_UNDEF -1.0
#define ARG_DEPTH_FILE -2.0
#define ARG_DEPTH_INF -3.0

#define ARG_ERROR_RATE_UNDEF -1.0

#define ARG_DEPTH_MAXVAL 500.0 

#define ARG_NTHREADS_MAXVAL 50

#define ARG_BETA_VAR_UNDEF -1.0


/// [args->doUnobserved]
/// notations for representing non-reference unobserved alleles

/// trim to include only observed bases
#define ARG_DOUNOBSERVED_TRIM 0

/// use <*> symbolic alternate allele (used in bcftools)
#define ARG_DOUNOBSERVED_STAR 1

/// use <NON_REF> gVCF NON_REF notation (used in GATK)
#define ARG_DOUNOBSERVED_NONREF 2

/// explode to unseen bases
#define ARG_DOUNOBSERVED_EXPLODE_ACGT 3

/// explode and add <*> to the end
#define ARG_DOUNOBSERVED_EXPLODE_ACGT_STAR 4 

/// explode and add <NON_REF> to the end
#define ARG_DOUNOBSERVED_EXPLODE_ACGT_NONREF 5 

/// retain REF/ALT alleles in the input file
#define ARG_RETAIN_REFALT 1

/// define the type of simulation source GT tag to use
/// REF     ALT    GT
/// 0       1      0/1
/// output from msprime BinaryMutationModel
#define ARG_GTSOURCE_BINARY 0

/// REF     ALT    GT
/// C       T,G      0/1
#define ARG_GTSOURCE_ACGT 1

#define ARG_I16_MAPQ_DEFAULT (20)

// args->adjustQs
// NONE: do not perform any adjustment
#define ARG_QS_ADJUST_DISABLED       (0)
// [+1] use adjusted quality scores for genotype likelihoods
#define ARG_QS_ADJUST_FOR_GL     (1<<0)
// [+2] use adjusted quality scores for qsum (QS tag)
#define ARG_QS_ADJUST_FOR_QSUM   (1<<1)
// [+4] use adjusted quality scores for pileup output
#define ARG_QS_ADJUST_FOR_PILEUP (1<<2)
// [+8] use adjusted quality scores for --printQScores
#define ARG_QS_ADJUST_FOR_PRINTQSCORES (1<<3)
// [+16] use adjusted quality scores for --printGlError
#define ARG_QS_ADJUST_FOR_PRINTGLERROR (1<<4)

#define ARGMAX_QS_ADJUST (ARG_QS_ADJUST_FOR_GL | ARG_QS_ADJUST_FOR_QSUM | ARG_QS_ADJUST_FOR_PILEUP | ARG_QS_ADJUST_FOR_PRINTQSCORES | ARG_QS_ADJUST_FOR_PRINTGLERROR)

#define ARG_RM_INVAR_DISABLED (0<<0)
#define ARG_RM_INVAR_INPUT_HOMOREFGT (1<<0)
#define ARG_RM_INVAR_INPUT_HOMOALTGT (1<<1)
#define ARG_RM_INVAR_INPUT_HOMOGT (ARG_RM_INVAR_INPUT_HOMOREFGT | ARG_RM_INVAR_INPUT_HOMOALTGT)
#define ARG_RM_INVAR_INPUT_SIM_INVAR (1<<2)

#define ARGMAX_RM_INVAR (ARG_RM_INVAR_INPUT_HOMOREFGT | ARG_RM_INVAR_INPUT_HOMOALTGT | ARG_RM_INVAR_INPUT_SIM_INVAR)


/// --> CONSTANTS:SIGNALS

/// int prepare_gvcf_block()
/// at bcf_utils.cpp
#define GVCF_NO_WRITE 0
#define GVCF_FLUSH_BLOCK 1
#define GVCF_WRITE_SIMREC 2



/* -> STATEMENTS --------------------------------------------------------------*/

/// --> STATEMENTS:PROGRAM WILL

#define PROGRAM_WILL(arg, flag) \
    ( ((arg) & (1<<(flag))) )

#define PROGRAM_BETA_VAR_IS_UNDEF \
    ( ((args->beta_variance) == (ARG_BETA_VAR_UNDEF)) )

#define PROGRAM_WILL_EXPLODE_ACGT \
    (((args->doUnobserved)==ARG_DOUNOBSERVED_EXPLODE_ACGT) || ((args->doUnobserved)==ARG_DOUNOBSERVED_EXPLODE_ACGT_STAR) || ((args->doUnobserved)==ARG_DOUNOBSERVED_EXPLODE_ACGT_NONREF))

#define PROGRAM_WILL_ADD_UNOBSERVED \
    ( (args->doUnobserved==ARG_DOUNOBSERVED_STAR) || (args->doUnobserved==ARG_DOUNOBSERVED_NONREF) || (args->doUnobserved==ARG_DOUNOBSERVED_EXPLODE_ACGT_STAR) || (args->doUnobserved==ARG_DOUNOBSERVED_EXPLODE_ACGT_NONREF) )

#define PROGRAM_WILL_ADD_STAR \
    ( (args->doUnobserved==ARG_DOUNOBSERVED_STAR) || (args->doUnobserved==ARG_DOUNOBSERVED_EXPLODE_ACGT_STAR) )

#define PROGRAM_WILL_ADD_NONREF \
    ( (args->doUnobserved==ARG_DOUNOBSERVED_NONREF) || (args->doUnobserved==ARG_DOUNOBSERVED_EXPLODE_ACGT_NONREF) )

#define PROGRAM_WILL_SAMPLE_STRAND \
    ( (args->addI16 || args->addFormatADF || args->addFormatADR || args->addInfoADF || args->addInfoADR))

#define PROGRAM_WILL_SKIP_INPUT_HOMOREFGT_SITES \
    ( ((args->rmInvarSites) & ARG_RM_INVAR_INPUT_HOMOREFGT) )

#define PROGRAM_WILL_SKIP_INPUT_HOMOALTGT_SITES \
    ( ((args->rmInvarSites) & ARG_RM_INVAR_INPUT_HOMOALTGT) )

#define PROGRAM_WILL_SKIP_INPUT_HOMOGT_SITES \
    ( ((args->rmInvarSites) & ARG_RM_INVAR_INPUT_HOMOGT) )

#define PROGRAM_WILL_SKIP_SIM_INVAR_SITES \
    ( ((args->rmInvarSites) & ARG_RM_INVAR_INPUT_SIM_INVAR) )

#define PROGRAM_WILL_USE_BINARY_GTSOURCE \
    ( ((args->gtSource)==ARG_GTSOURCE_BINARY) )

#define PROGRAM_WILL_USE_ACGT_GTSOURCE \
    ( ((args->gtSource)==ARG_GTSOURCE_ACGT) )

#define PROGRAM_WILL_ADJUST_QS_FOR_GL \
    ( ((args->adjustQs) & ARG_QS_ADJUST_FOR_GL) )

#define PROGRAM_WILL_ADJUST_QS_FOR_QSUM \
    ( ((args->adjustQs) & ARG_QS_ADJUST_FOR_QSUM) )

#define PROGRAM_WILL_ADJUST_QS_FOR_PILEUP \
    ( ((args->adjustQs) & ARG_QS_ADJUST_FOR_PILEUP) )

#define PROGRAM_WILL_ADJUST_QS_FOR_PRINTQSCORES \
    ( ((args->adjustQs) & ARG_QS_ADJUST_FOR_PRINTQSCORES) )

#define PROGRAM_WILL_ADJUST_QS_FOR_PRINTGLERROR \
    ( ((args->adjustQs) & ARG_QS_ADJUST_FOR_PRINTGLERROR) )

#define PROGRAM_WILL_ADJUST_QS \
    ( (args->adjustQs) )

#define BASE_A 0
#define BASE_C 1
#define BASE_G 2
#define BASE_T 3
#define BASE_NONREF 4

#define MAXGL 0.0 /// best
#define MINGL NEG_INF /// worst 

#define MAXPL 255 /// worst
#define MINPL 0 /// best

#define MAXGP 1.0 /// best
#define MINGP 0.0 /// worst

#define SIM_FORWARD_STRAND 0
#define SIM_REVERSE_STRAND 1

/// maximum number of alleles to consider 
/// | 0 | 1 | 2 | 3 | 4 |
/// | A | C | G | T |<*>|
#define MAX_NALLELES 5  

/// maximum number of genotypes to consider
/// depends on: MAX_NALLELES
///  (MAX_NALLELES * (MAX_NALLELES+1)) / 2
/// |  0  |  1  |  5  |  2  |  6  |  9  |  3  |  7  |  10 |  12 | 4      | 8     | 11    | 13    | 14  		|
/// |  AA |  AC |  CC |  AG |  CG |  GG |  AT |  CT |  GT |  TT | A<*> 	| C<*> 	| G<*> 	| T<*> 	| <*><*> 	|
/// A0A0 A0A1 A1A1 A0A2 A1A2 A2A2 A0A3 A1A3 A2A3 A3A3 A0A4 A1A4 A2A4 A3A4 A4A4
#define MAX_NGTS 15  

/// Old order (vcfgl v<=0.3.3)
/// AA, AC, AG, AT, A<*>, CC, CG, CT, C<*>, GG, GT, G<*>, TT, T<*>, <*><*>
/// AA, AC, AG, AT, CC, CG, CT, GG, GT, TT

/// assume diploid samples
#define SIM_PLOIDY 2

/// source: bcftools/bam2bcf.c L41 CAP_DIST
#define CAP_TAIL_DIST 25

/// source: bcftools/bam2bcf.c L381
#define CAP_BASEQ 63

#define BCF_GT_PHASED_0 3  /// bcf_gt_phased(0)
#define BCF_GT_PHASED_1 5  /// bcf_gt_phased(1)


/* -> FUNCTION-LIKE MACROS ---------------------------------------------------*/

#define LOG2LN(x) ((((x)) / (M_LOG10E)))


#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))


#define FLUSH_BGZF_KSTRING_BUFFER(fp, buffer) \
    do {                                               \
        if ( (((buffer)->l) > KSTRING_BGZF_WRITE_BUFFER_SIZE) ) { \
            write_BGZF(fp, (buffer)->s, (buffer)->l);  \
            (buffer)->l = 0;                           \
        } \
    } while (0);

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
#define ASSERT_EXPAND(x) x
#define ASSERT(expr)                                                        \
    do {                                                                    \
        if (!(ASSERT_EXPAND(expr))){ \
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


#define CHECK_ARG_INTERVAL_INT(argval, minval, maxval, argstr) \
do { \
    if (((argval) < (minval)) || ((argval) > (maxval)) ) { \
        \
            ERROR("[Bad argument value: '%s %d'] Allowed range is [%d,%d]", (argstr), (argval), (minval), (maxval)); \
    } \
} while (0);


#define CHECK_ARG_INTERVAL_DBL(argval, minval, maxval, argstr) \
do { \
    if (((argval) < (minval)) || ((argval) > (maxval)) ) { \
        \
            ERROR("[Bad argument value: '%s %f'] Allowed range is [%f,%f]", (argstr), (argval), (minval), (maxval)); \
    } \
} while (0);

#define CHECK_ARG_INTERVAL_IE_DBL(argval, minval, maxval, argstr) \
do { \
    if (((argval) < (minval)) || ((argval) >= (maxval)) ) { \
        \
            ERROR("[Bad argument value: '%s %f'] Allowed range is [%f,%f]", (argstr), (argval), (minval), (maxval)); \
    } \
} while (0);


#define CHECK_ARG_INTERVAL_II_DBL(argval, minval, maxval, argstr) \
do { \
    if (((argval) < (minval)) || ((argval) > (maxval)) ) { \
        \
            ERROR("[Bad argument value: '%s %f'] Allowed range is [%f,%f]", (argstr), (argval), (minval), (maxval)); \
    } \
} while (0);


#define CHECK_ARG_INTERVAL_EI_DBL(argval, minval, maxval, argstr) \
do { \
    if (((argval) <= (minval)) || ((argval) > (maxval)) ) { \
        \
            ERROR("[Bad argument value: '%s %f'] Allowed range is [%f,%f]", (argstr), (argval), (minval), (maxval)); \
    } \
} while (0);


#define CHECK_ARG_INTERVAL_EE_DBL(argval, minval, maxval, argstr) \
do { \
    if (((argval) <= (minval)) || ((argval) >= (maxval)) ) { \
        \
            ERROR("[Bad argument value: '%s %f'] Allowed range is [%f,%f]", (argstr), (argval), (minval), (maxval)); \
    } \
} while (0);


#define CHECK_ARG_INTERVAL_01(argval, argstr) \
do { \
    if ( ((argval)!=0) && ((argval)!=1) ) { \
        \
            ERROR("Argument %s with value %d is out of range. Allowed values are 0 (for on/enable) and 1 (for off/disable)", (argstr), (argval)); \
    } \
} while (0);

#define CHECK_ARG_VALUES_LIST(argval, argstr, ...) \
do { \
    int __argval = (argval); \
    int __argvals[] = { __VA_ARGS__ }; \
    int __nvals = sizeof(__argvals) / sizeof(int); \
    int __i; \
    int __found = 0; \
    for (__i = 0; __i < __nvals; ++__i) { \
        if (__argval == __argvals[__i]) { \
            __found = 1; \
            break; \
        } \
    } \
    if (!__found) { \
        ERROR("Bad argument value: '%s %d'", (argstr), (argval)); \
    } \
} while (0);



           /* LOOKUP TABLES & LOOKUP MACROS ------------------------------------------*/

extern const int lut_char_to_acgt_int[85];

#define ACGT_CHAR2INT(c) \
    ( ( (((unsigned char)(c)) < 85) ? (lut_char_to_acgt_int[((unsigned char)(c))]) : (-1) ) )


/// @brief lut_qs_to_qs2 - map quality score to squared quality score
/// @details
/// [R]
/// > n <- 63
/// > paste0("const int qs_to_qs2[", n + 1, "] = {", paste(sapply(0:n,
/// function(n) n*n),collapse = ", "),"};") [1] "const int qs_to_qs2[64] = {0, 1,
/// 4, 9, 16, 25, 36, 49, 64, 81, 100, 121, 144, 169, 196, 225, 256, 289, 324,
/// 361, 400, 441, 484, 529, 576, 625, 676, 729, 784, 841, 900, 961, 1024, 1089,
/// 1156, 1225, 1296, 1369, 1444, 1521, 1600, 1681, 1764, 1849, 1936, 2025, 2116,
/// 2209, 2304, 2401, 2500, 2601, 2704, 2809, 2916, 3025, 3136, 3249, 3364, 3481,
/// 3600, 3721, 3844, 3969};"
/// remove CAP_BASEQ (63)
extern const int lut_qs_to_qs2[63];
/// @assume CAP_BASEQ == 63
#define QS_TO_QSSQ(q) ( (0==(q)) ? (0) : ( ((q)<CAP_BASEQ) ? (lut_qs_to_qs2[(q)]) : (3969) ) )


/// @brief NALLELES_TO_NGTS - get number of genotypes from number of
/// alleles
/// @param 	n (int) - number of alleles
/// @return	  (int) - number of unique unordered genotypes expected for n
/// alleles
///                    e.g. A, C -> AA AC CC -> 3
extern const int lut_nAlleles_to_nGenotypes[6];
#define NALLELES_TO_NGTS(n) ( (((n)>0) && ((n)<6)) ? lut_nAlleles_to_nGenotypes[(n)] : -1 )

/// @note not in use - alternative method
/// @brief nAlleles2nGenotypes - get the number of possible genotypes assuming
/// ploidy==2
/// @param n	
/// number of alleles
/// @return		number of expected genotypes
/// equivalent to (n * (n+1)) / 2
inline int nAlleles2nGenotypes(const int n) {
    return (((n * (n + 1)) >> 1));
}



/// [R]
/// qToP<-function(q){10^(-q/10)}
/// CAP_BASEQ=63
/// arrSize=CAP_BASEQ+1
/// paste0("qScore_to_errorProb[",arrSize,"]={",paste(unlist(format(lapply(0:CAP_BASEQ,FUN=qToP),scientific=F)),collapse=","),"};")
/// remove CAP_BASEQ (63)
extern const double qScore_to_errorProb[63];
/// @assume CAP_BASEQ == 63
/// @note does NOT check if q is in range [0, CAP_BASEQ]; DRAGON
#define QS_TO_ERRPROB(q) ( (0==(q)) ? 1.0 : ( ((q) < CAP_BASEQ) ? qScore_to_errorProb[(q)] : 0.0000005011872 ) )

/// @brief qScore_to_ln_gl - lookup table for mapping quality score to ln genotype likelihoods
/// [R]
/// qToP <- function(q) { 10 ^ (-q / 10) }
/// CAP_BASEQ = 256
/// arrSize = CAP_BASEQ + 1
/// generate_likelihoods_homhit_ln <- function(i) { p = 1 - qToP(i);ret = log((p * 0.5) + (p * 0.5));return(ret); }
/// generate_likelihoods_hethit_ln <- function(i) { p = 1 - qToP(i);p3 = (1 - p) / 3.0;ret = log((p * 0.5) + (p3 * 0.5));return(ret); }
/// generate_likelihoods_homnonhit_ln <- function(i) { p = 1 - qToP(i);p3 = (1 - p) / 3.0;ret = log((0.5 * p3) + (0.5 * p3));return(ret); }
/// fn = generate_likelihoods_homhit_ln
/// x1 = paste0("{", paste(unlist(format(lapply(0:CAP_BASEQ, FUN = fn), scientific = F)), collapse = ","), "}")
/// fn = generate_likelihoods_hethit_ln
/// x2 = paste0("{", paste(unlist(format(lapply(0:CAP_BASEQ, FUN = fn), scientific = F)), collapse = ","), "}")
/// fn = generate_likelihoods_homnonhit_ln
/// x3 = paste0("{", paste(unlist(format(lapply(0:CAP_BASEQ, FUN = fn), scientific = F)), collapse = ","), "}")
/// print(gsub("-Inf", "NEG_INF", paste0("qScore_to_ln_gl[3][", arrSize, "]=", "{", paste(x1, x2, x3, sep = ","), "};")))
extern const double qScore_to_ln_gl[3][257];

/// @brief qScore_to_log10_gl - lookup table for mapping quality score to log10 genotype likelihoods
/// [R]
/// qToP <- function(q) { 10 ^ (-q / 10) }
/// CAP_BASEQ = 256
/// arrSize = CAP_BASEQ + 1
/// generate_likelihoods_homhit_log10 <- function(i) { p = qToP(i);ret = log10(1 - p);return(ret); }
/// generate_likelihoods_hethit_log10 <- function(i) { p = qToP(i);ret = log10((1 - p) / 2 + p / 6);return(ret); }
/// generate_likelihoods_homnonhit_log10 <- function(i) { p = qToP(i);ret = log10(p) - log10(3);return(ret); }
/// fn = generate_likelihoods_homhit_log10
/// x1 = paste0("{", paste(unlist(format(lapply(0:CAP_BASEQ, FUN = fn), scientific = F)), collapse = ","), "}")
/// fn = generate_likelihoods_hethit_log10
/// x2 = paste0("{", paste(unlist(format(lapply(0:CAP_BASEQ, FUN = fn), scientific = F)), collapse = ","), "}")
/// fn = generate_likelihoods_homnonhit_log10
/// x3 = paste0("{", paste(unlist(format(lapply(0:CAP_BASEQ, FUN = fn), scientific = F)), collapse = ","), "}")
/// print(gsub("-Inf", "NEG_INF", paste0("qScore_to_log10_gl[3][", arrSize, "]=", "{", paste(x1, x2, x3, sep = ","), "};")))
extern const double qScore_to_log10_gl[3][257];

/// @brief qScore_to_gl - lookup table for mapping quality score to genotype likelihoods
/// [R]
/// qToP <- function(q) { 10 ^ (-q / 10) }
/// CAP_BASEQ = 256
/// arrSize = CAP_BASEQ + 1
/// generate_likelihoods_homhit_ln < -function(i) { p = 1 - qToP(i);ret = ((p * 0.5) + (p * 0.5));return(ret); }
/// generate_likelihoods_hethit_ln < -function(i) { p = 1 - qToP(i);p3 = (1 - p) / 3.0;ret = ((p * 0.5) + (p3 * 0.5));return(ret); }
/// generate_likelihoods_homnonhit_ln < -function(i) { p = 1 - qToP(i);p3 = (1 - p) / 3.0;ret = ((0.5 * p3) + (0.5 * p3));return(ret); }
/// fn = generate_likelihoods_homhit_ln
/// x1 = paste0("{", paste(unlist(format(lapply(0:CAP_BASEQ, FUN = fn), scientific = F)), collapse = ","), "}")
/// fn = generate_likelihoods_hethit_ln
/// x2 = paste0("{", paste(unlist(format(lapply(0:CAP_BASEQ, FUN = fn), scientific = F)), collapse = ","), "}")
/// fn = generate_likelihoods_homnonhit_ln
/// x3 = paste0("{", paste(unlist(format(lapply(0:CAP_BASEQ, FUN = fn), scientific = F)), collapse = ","), "}")
/// print(gsub("-Inf", "NEG_INF", paste0("qScore_to_gl[3][", arrSize, "]=", "{", paste(x1, x2, x3, sep = ","), "};")))
extern const double qScore_to_gl[3][257];

extern float bcf_float_missing_union_f;



#endif  /// __SHARED__
