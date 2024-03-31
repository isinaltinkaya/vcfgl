#ifndef __ARGUMENTS__
#define __ARGUMENTS_

#include <htslib/kstring.h>  // kstring_t
#include <htslib/bgzf.h> // bgzf
#include <htslib/hts.h>  // hts_version()

#include <ctime> // time_t, localtime, asctime

#include "random_generator.h" // BetaSampler, PoissonSampler


FILE* get_FILE(const char* fname, const char* mode);
FILE* open_FILE(const char* a, const char* b);
htsFile* open_htsFile(const char* fn, const char* mode);
BGZF* open_BGZF(const char* fn, const char* mode);
void write_BGZF(BGZF* fp, const void* data, const int size);
void write_BGZF_kstring_buffer(BGZF* fp, kstring_t* buffer);


typedef struct preCalcStruct preCalcStruct;

struct preCalcStruct {
    double homT = -1.0;
    double het = -1.0;
    double homF = -1.0;
    double error_prob_forQs = -1.0;
    double error_prob_forGl = -1.0;
    int qScore = -1;
    int q5 = -1; // qScore << 5 value
    int adj_qScore = -1;
    int adj_q5 = -1; // adj_qScore << 5 value
};


extern void help_page();

typedef struct argStruct argStruct;
extern argStruct* args;

struct argStruct {

    // -------------------------------------------- //
    // read in from command line

    int verbose;
    int n_threads;
    int seed;

    // --> Input/Output

    char* in_fn;
    int gtSource;
    char* out_fnprefix;
    char* output_mode;

    // --> Simulation parameters

    // mps_depth: mean per-site read depth
    // ARG_DEPTH_ISINF     (--depth "inf")
    // ARG_DEPTH_ISFILE       (--depths-file <file>)
    // x        (--depth x)
    // ARG_DEPTH_ISUNDEF (default)
    double mps_depth;
    // mps_depths_fn: file containing mean per-site read depth for each sample
    char* mps_depths_fn;
    double error_rate;
    int error_qs;
    double beta_variance;
    int GL;
    double glModel1_theta; // err dep parameter
    char* qs_bins_fn;
    int usePreciseGlError;
    int i16_mapq;
    char* gvcf_dps_str;
    int adjustQs;
    double adjustBy;

    int explode;
    int rmInvarSites;
    int rmEmptySites;
    int doUnobserved;
    int doGVCF;
    int printPileup;
    int printTruth;
    int printBasePickError;
    int printQsError;
    int printGlError;
    int printQScores;

    int addGL;
    int addGP;
    int addPL;
    int addI16;
    int addQS;
    int addFormatDP;
    int addInfoDP;
    int addFormatAD;
    int addInfoAD;
    int addFormatADF;
    int addInfoADF;
    int addFormatADR;
    int addInfoADR;

    // -------------------------------------------- //
    // set during program execution

    char* out_fn;
    char* out_truth_fn;
    char* out_pileup_fn;
    char* output_mode_str;

    htsFile* in_fp;
    FILE* arg_fp;
    htsFile* out_fp;
    htsFile* out_truth_fp;
    BGZF* out_pileup_fp;

    char* datetime;
    char* command;
    char* versionInfo;

    // -------------------------------------------- //
    // for internal use

    double* mps_depths;

    int n_mps_depths; // 1 if mps_depths_fn == NULL, otherwise n_mps_depths = n_samples
    BetaSampler* betaSampler;
    PoissonSampler** poissonSampler;

    double base_pick_error_prob;

    // precalculation mode on/off
    // if error_qs == 0, then preCalc !=NULL
    // if error_qs == 1, then preCalc !=NULL
    // if error_qs == 2, then preCalc ==NULL
    preCalcStruct* preCalc;

    errmod_t* gl1errmod;

    /// qs_bins[nRanges]
    /// qs_bins[i][0] = start of i-th range
    /// qs_bins[i][1] = end of i-th range
    /// qs_bins[i][2] = qs value to assign for i-th range
    uint8_t** qs_bins;
    int n_qs_bins;

    // functions:
    double* read_depthsFile(void);


};

argStruct* args_init();

argStruct* args_get(int argc, char** argv);

void args_destroy(argStruct* args);


#endif  // __ARGUMENTS__
