#ifndef __ARGUMENTS__
#define __ARGUMENTS_

#include <htslib/kstring.h> // kstring_t
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <sys/stat.h>
#include <stdio.h>

#include "version.h"

#include "shared.h"

#include <htslib/hts.h> // hts_version()

#include "random_generator.h"

FILE *getFILE(const char *fname, const char *mode);
FILE *openFILE(const char *a, const char *b);

kstring_t *kbuf_init();
void kbuf_destroy(kstring_t *kbuf);

extern void help_page();

typedef struct argStruct argStruct;
extern argStruct *args;

/*
 *
 * @field *datetime			date and time of execution
 * @field *command			command line
 *
 * @field *in_fn			input filename
 * @field *out_fnp			output filename prefix
 * @field error_rate			error rate
 * @field mps_depth			mean per-site read depth
 * @field mps_depths_fn		assign depths to individuals from per individual mean per site depth file, one line per individual
 *
 * @field error_bias		should the program sample errors?
 * 							0: no
 * 							1: sample error probabilities from beta distribution with mean=error_rate var=beta_variance
 *
 * @field beta_variance	variance of the beta distribution
 *
 * @field pos0				are positions 0-based? [0]
 *							if 1 is set; input VCF positions are 0-based;
 *							shift coordinate system+1;
 * @field seed
 * @field output_mode		char defining the output file format
 *
 *							Output modes
 *						b	compressed BCF [default]
 *						u	uncompressed BCF
 *						v	uncompressed  VCF
 *						z	compressed VCF (bgzf compressed)
 *
 * @field explode			explode to unobserved invariable sites
 * @field printBaseCounts	should the program print base counts
 *
 * @field addGP				add GP field
 * @field addPL				add PL field
 * @field addI16			add I16 field
 * @field addQS				add QS field
 *
 *
 */
struct argStruct
{

	char **argv;

	char *datetime;
	char *command;

	BetaSampler *betaSampler = NULL;

	char *in_fn;
	char *out_fnp;

	double error_rate;
	double mps_depth;
	char *mps_depths_fn;

	int error_bias;
	double beta_variance;
	int error_rate_q;

	int pos0;
	int seed;
	int trimAlts;

	int addGP;
	int addPL;
	int addI16;
	int addQS;

	char *output_mode;

	int explode;
	int printBaseCounts;
};

argStruct *args_init();

argStruct *args_get(int argc, char **argv);

double *read_depthsFile(const char *fname, int len);

void args_destroy(argStruct *args);

#endif // __ARGUMENTS__
