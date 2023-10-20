#ifndef __ARGUMENTS__
#define __ARGUMENTS_

#include <ctype.h>
#include <htslib/kstring.h>  // kstring_t
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

#include "version.h"

#include "shared.h"

#include <htslib/hts.h>  // hts_version()

#include "random_generator.h"

FILE* getFILE(const char* fname, const char* mode);
FILE* openFILE(const char* a, const char* b);

extern void help_page();

typedef struct argStruct argStruct;
extern argStruct* args;

struct argStruct {

    // -------------------------
    // read in from command line

    int verbose;
    int n_threads;

    char* in_fn;
    char* out_fnp;
    char* output_mode;

    // mps_depth: mean per-site read depth
    // -999.0     (--depth "inf")
    // -1.0       (--depths-file <file>)
    // x        (--depth x)
    double mps_depth;
    char* mps_depths_fn;

    double error_rate;
    int error_qs;
    double beta_variance;

    int usePreciseGlError;


    int pos0;
    int seed;

    int trimAlts;
    int rmInvarSites;
    int useUnknownAllele;

    int platform;
    int explode;

    int addGL;
    int addGP;
    int addPL;
    int addI16;
    int addQS;
    int addFormatDP;
    int addFormatAD;

    int addFormatADF;
    int addFormatADR;

    int addInfoDP;
    int addInfoAD;

    int addInfoADF;
    int addInfoADR;

    // ----------------------------
    // set during program execution

    char* out_fn;
    char* output_mode_str;

    char* datetime;
    char* command;

    BetaSampler* betaSampler = NULL;

    double base_pick_error_prob;
    double error_prob_forQs; // e for qs calculation
    int qScore;
    double error_prob_forGl; // e for gl calculation; interacts with usePreciseGlError

    double homT;
    double het;
    double homF;

    FILE* arg_ff;

};

argStruct* args_init();

argStruct* args_get(int argc, char** argv);

double* read_depthsFile(const char* fname, int len);

void args_destroy(argStruct* args);

#endif  // __ARGUMENTS__
