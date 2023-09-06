#ifndef __SHARED__
#define __SHARED__


#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <limits>
#include "dev.h"

/* ========================================================================== */
/* /BEGIN/ MACRO DEFINITIONS =================================================*/
/* ========================================================================== */

/* -> CONSTANTS --------------------------------------------------------------*/

#define MAXGL 0
#define MINGL NEG_INF
#define MAXPL 255
#define MINPL 0
#define MAXGP 1.0
#define MINGP 0.0

// number of genotypes to simulate {AA,AC,CC,AG,CG,GG,AT,CT,GT,TT}
#define SIM_NGTS 10

#define SIM_PLOIDY 2

// source: bcftools/bam2bcf.c L41
#define CAP_DIST 25

// source: bcftools/bam2bcf.c L381
#define CAP_BASEQ 63

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
#define ERROR(...)                                                                           \
	do{ \
		fprintf(stderr, "\n\n*******\n[ERROR](%s/%s:%d)\n\t", __FILE__, __FUNCTION__, __LINE__); \
		fprintf(stderr, __VA_ARGS__);                                                            \
		fprintf(stderr, "\n*******\n");                                                          \
		exit(1); \
	}while(0);

/*
 * Macro:[NEVER]
 * indicates that a point in the code should never be reached
 */
#define NEVER \
	do{ \
		ERROR("Control should never reach this point; please report this to the developers.") \
	}while(0);

/*
 * Macro:[ASSERT]
 * evaluate an expression, works the same way as the C-macro assert
 * except that DEBUG does not affect it (it is always active)
 * also prints the file and line info and exits the program
 * if the expression evaluates to false
 */
#define ASSERT(expr)                                                                                              \
	do{ \
		if (!((expr))) {                                                                                              \
			fprintf(stderr, "\n\n*******\n[ERROR](%s/%s:%d) %s\n*******\n", __FILE__, __FUNCTION__, __LINE__, #expr); \
			exit(1);                                                                                                  \
		} \
	}while(0);


/*
 * Macro:[WARNING]
 * print a custom warning message
 */
#define WARNING(...)                                                                \
	do{ \
		fprintf(stderr, "\n\n[WARNING](%s/%s:%d): ", __FILE__, __FUNCTION__, __LINE__); \
		fprintf(stderr, __VA_ARGS__);                                                   \
		fprintf(stderr, "\n");                                                          \
	}while(0);


/* ========================================================================== */
/* /END/ MACRO DEFINITIONS ================================================== */
/* ========================================================================== */

const double NEG_INF = -std::numeric_limits<double>::infinity();






#endif // __SHARED__
