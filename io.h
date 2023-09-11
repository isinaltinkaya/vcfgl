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

FILE *getFILE(const char *fname, const char *mode);
FILE *openFILE(const char *a, const char *b);

kstring_t *kbuf_init();
void kbuf_destroy(kstring_t *kbuf);

extern void help_page();

typedef struct argStruct argStruct;
extern argStruct *args;

/*
 *
 * @field *in_fn			pointer to input file name
 * @field *out_fp			pointer to output file prefix
 * @field errate			error rate [0.01]
 * @field mps_depth			mean per site depth [1]
 * @field in_mps_depths		assign depths to individuals from per individual mean per site depth file, one line per individual
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
 * @field explode
 * @field printBaseCounts	should the program print base counts
 */
struct argStruct
{

	char **argv;

	char *in_fn;
	char *out_fp;

	char *datetime;
	char *command;

	double errate;
	double mps_depth;
	char *in_mps_depths;

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
