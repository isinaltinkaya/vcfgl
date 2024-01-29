#include "io.h"
#include "gl_methods.h"
#include "version.h"
#include "build.h"
#include <time.h> // asctime
#include <sys/stat.h> // stat()


extern void(*calculate_gls)(simRecord* sim);
extern unsigned short int rng1_seeder[3];
extern unsigned short int rng1_seeder_save[3];
extern unsigned short int rng2_seeder[3];
extern unsigned short int rng2_seeder_save[3];
extern const char* nonref_str;


preCalcStruct::preCalcStruct() {
    this->homT = -1.0;
    this->het = -1.0;
    this->homF = -1.0;
    this->error_prob_forQs = -1.0;
    this->error_prob_forGl = -1.0;
    this->qScore = -1;
    this->q5 = -1;
}

void preCalcStruct::prepare_gls_preCalc(void) {

    if (0 == args->usePreciseGlError) {

        this->error_prob_forGl = QS_TO_ERRPROB(this->qScore);
        this->homT = qScore_to_log10_gl[0][qScore];
        this->het = qScore_to_log10_gl[1][qScore];
        this->homF = qScore_to_log10_gl[2][qScore];

    } else if (1 == args->usePreciseGlError) {

        this->error_prob_forGl = this->error_prob_forQs;
        if (0.0 == this->error_prob_forGl) {
            this->homT = 0;
            this->het = -0.30103;
            this->homF = NEG_INF;
        } else if (this->error_prob_forGl > 0.0) {
            // prepare for gl calculation with pre-calculated terms
            this->homT = log10(1.0 - this->error_prob_forGl);
            this->het = log10((1.0 - this->error_prob_forGl) / 2.0 + this->error_prob_forGl / 6.0);
            this->homF = log10(this->error_prob_forGl) - PRE_CALC_LOG10_3;
        } else {
            NEVER;
        }
    }
}


FILE* get_FILE(const char* fname, const char* mode) {
    FILE* fp;
    if (NULL == (fp = fopen(fname, mode))) {
        fprintf(stderr,
            "[%s:%s()]\t->Error opening FILE handle for file:%s exiting\n",
            __FILE__, __FUNCTION__, fname);
        exit(0);
    }
    return fp;
}

FILE* open_FILE(const char* prefix, const char* suffix) {
    char* str = (char*)malloc(strlen(prefix) + strlen(suffix) + 1);
    strcpy(str, prefix);
    strcat(str, suffix);
    // fprintf(stderr, "\t-> Opening output file: %s\n", str);
    FILE* fp = get_FILE(str, "w");
    free(str);
    str = NULL;
    return(fp);
}

htsFile* open_htsFile(const char* fn, const char* mode) {
    htsFile* fp = hts_open(fn, mode);
    if (NULL == fp) {
        ERROR("Could not open file: %s with mode %s", fn, mode);
    }
    return(fp);
}

BGZF* open_BGZF(const char* fn, const char* mode) {
    BGZF* fp = bgzf_open(fn, mode);
    if (NULL == fp) {
        ERROR("Could not open file: %s with mode %s", fn, mode);
    }
    return(fp);
}

void write_BGZF(BGZF* fp, const void* data, const int size) {
    if ((bgzf_write(fp, data, size)) != size) {
        ERROR("Could not write to BGZF file");
    }
    return;
}


char* get_time() {
    time_t current_time;
    struct tm* local_time;
    current_time = time(NULL);
    local_time = localtime(&current_time);
    return (asctime(local_time));
}


void version_page() {
    fprintf(stderr, "vcfgl [version: %s] [build: %s %s] [htslib: %s]\n", VCFGL_VERSION, __DATE__, __TIME__, hts_version());
    fprintf(stderr, "\n");
    fprintf(stderr, "Build details:\n");
    fprintf(stderr, "\t-> CXX=%s\n", VCFGL_MAKE_CXX);
    fprintf(stderr, "\t-> CXXFLAGS=%s\n", VCFGL_MAKE_CXXFLAGS);
    fprintf(stderr, "\t-> CPPFLAGS=%s\n", VCFGL_MAKE_CPPFLAGS);
    fprintf(stderr, "\t-> LIBS=%s\n", VCFGL_MAKE_LIBS);
    fprintf(stderr, "\t-> FLAGS=%s\n", VCFGL_MAKE_FLAGS);
    fprintf(stderr, "\t-> HTSSRC=%s\n", VCFGL_MAKE_HTSSRC);
    fprintf(stderr, "\n");

}

void help_page() {
    // fprintf(stderr, "vcfgl [version: %s] [build: %s %s] [htslib: %s]\n", VCFGL_VERSION, __DATE__, __TIME__, hts_version());

    fprintf(stderr, "\n");
    fprintf(stderr, "Program: vcfgl\n");
    fprintf(stderr, "License: GNU GPLv3.0\n");
    fprintf(stderr, "Version: %s (htslib: %s)\n", VCFGL_VERSION, hts_version());
    fprintf(stderr, "Build: %s %s\n", __DATE__, __TIME__);
    fprintf(stderr, "\n");

    fprintf(stderr, "Usage: vcfgl -i <input> [options]\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "    -h, --help _____________________  Print this help message and exit\n");
    fprintf(stderr, "    -v, --version __________________  Print version and build information and exit\n");
    fprintf(stderr, "\n");

    fprintf(stderr, "Option descriptions:\n");
    fprintf(stderr, "     -s, --long-option TYPE [X] _____ Description\n");
    fprintf(stderr, "     -s                               Short option (if any)\n");
    fprintf(stderr, "         --long-option                Long option\n");

    fprintf(stderr, "                       TYPE           Type of the argument value, can be:\n");
    fprintf(stderr, "                                        - INT (integer)\n");
    fprintf(stderr, "                                        - INT+ (additive integer: sum values to use together\n");
    fprintf(stderr, "                                        - FLOAT (floating point number)\n");
    fprintf(stderr, "                                        - STRING (string)\n");
    fprintf(stderr, "                                        - FILE (filename)\n");
    fprintf(stderr, "                                        - x|y|z (one of the listed values x, y or z)\n");
    fprintf(stderr, "                            [X]       Default argument value (if any)\n");
    fprintf(stderr, "                                _____ Connector to the option description for better readability\n");
    fprintf(stderr, "\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "General options:\n");
    fprintf(stderr, "    -V, --verbose INT [0] ___________ Verbosity level\n");
    fprintf(stderr, "    -@, --threads INT [1] ___________ Number of threads\n");
    fprintf(stderr, "    -s, --seed INT [time] ___________ Random seed for initializing the random number generator\n");
    fprintf(stderr, "\n");

    fprintf(stderr, "Input/Output:\n");
    fprintf(stderr, "    -i, --input FILE ________________ Input VCF/BCF file\n");
    fprintf(stderr, "        --source [0]|1 ______________ 0: Input REF/ALT alleles are in binary format (REF=0, ALT=1; typically outputted from msprime BinaryMutationModel)\n");
    fprintf(stderr, "                                      1: Input REF/ALT alleles are in VCF format (REF=i, ALT=j(,k..); i, j and k from {A,C,G,T}; i.e. the regular VCF format)\n");
    fprintf(stderr, "    -o, --output STRING ['output'] __ Output filename prefix\n");
    fprintf(stderr, "    -O, --output-mode [b]|u|z|v _____ b: Compressed BCF (.bcf), u: uncompressed BCF (.bcf), z: compressed VCF (.vcf.gz), v: uncompressed VCF (.vcf)\n");
    fprintf(stderr, "\n");

    fprintf(stderr, "Simulation parameters:\n");
    fprintf(stderr, "    -d, --depth FLOAT|'inf' _________ Mean per-site read depth\n");
    fprintf(stderr, "                                      ('inf') Simulate true values (requires: -addFormatDP 0 -addInfoDP 0)\n");
    fprintf(stderr, "   -df, --depths-file FILE __________ File containing mean per-site read depth values for each sample. One value per line.\n");
    fprintf(stderr, "    -e, --error-rate FLOAT __________ Base-calling error probability\n");
    fprintf(stderr, "   -eq, --error-qs [0]|1|2 __________ 0: Do not simulate errors in quality scores. Assumes all quality score assignments are correct\n");
    fprintf(stderr, "                                      1: Simulate site-specific errors in the probability of wrong base calls (requires: -bv FLOAT)\n");
    fprintf(stderr, "                                      2: Simulate the errors in the reported quality scores and genotype likelihoods (requires: -bv FLOAT)\n");
    fprintf(stderr, "   -bv, --beta-variance FLOAT _______ Designated variance for the beta distribution\n");
    fprintf(stderr, "   -GL, --gl-model 1|[2] ____________ Genotype likelihood model to be used in simulation\n");
    fprintf(stderr, "                                      1: Genotype likelihood model with correlated errors (a.k.a. samtools model, angsd -GL 1)\n");
    fprintf(stderr, "                                      2: Canonical genotype likelihood model with independent errors (a.k.a. GATK model, angsd -GL 2)\n");
    fprintf(stderr, "         --gl1-theta FLOAT [0.83] ___ Theta parameter for the genotype likelihood model 1 (requires: -GL 1)\n");
    fprintf(stderr, "         --platform [0]|1 ___________ 0: Do not use platform specification\n");
    fprintf(stderr, "                                      1: NovaSeq 6000 (qualities are binned into 4 values: 2, 12, 23 and 37)\n");

    fprintf(stderr, "         --precise-gl [0]|1 _________ 0: Use the discrete phred-scaled error probabilities in the genotype likelihood calculation\n");
    fprintf(stderr, "                                      1: Use precise error probabilities in the genotype likelihood calculation (requires: -GL 2)\n");
    fprintf(stderr, "         --i16-mapq INT [20] ________ Mapping quality score for I16 tag (requires: -addI16 1)\n");
    fprintf(stderr, "         --gvcf-dps INT(,INT..) _____ Minimum per-sample read depth range(s) for constructing gVCF blocks (requires: -doGVCF 1)\n");
    fprintf(stderr, "                                      Example: `--gvcf-dps 5,10,20` will group invariable sites into three types of gVCF blocks: [5,10), [10,20) and [20,inf).\n");
    fprintf(stderr, "                                      Sites with minimum depth < 5 will be printed as regular VCF records.\n");
    fprintf(stderr, "         -explode [0]|1 _____________ 1: Explode to sites that are not in input file.\n");
    fprintf(stderr, "                                      Useful for simulating invariable sites when the input file only contains variable sites.\n");
    fprintf(stderr, "                                      Sets all genotypes in exploded sites to homozygous reference.\n");
    fprintf(stderr, "         --rm-invar-sites INT+ ______ 0: Do not remove invariable sites\n");
    fprintf(stderr, "                                      1: Remove sites where all individuals' true genotypes in the input file are homozygous reference\n");
    fprintf(stderr, "                                      2: Remove sites where all individuals' true genotypes in the input file are homozygous alternative\n");
    fprintf(stderr, "                                      4: Remove sites where the all simulated reads among all individuals are the same base\n");
    fprintf(stderr, "                                      Example: '--rm-invar-sites 3' (1+2) will do both 1 and 2 (i.e. remove all homozygous sites)\n");
    fprintf(stderr, "         --rm-empty-sites [0]|1 _____ 0: Do not remove empty sites\n");
    fprintf(stderr, "                                      1: Remove empty sites (i.e. sites where no reads were simulated)\n");
    fprintf(stderr, "         -doUnobserved INT [1] ______ 0: Trim unobserved alleles. Only alleles that are observed will be listed in REF and ALT fields\n");
    fprintf(stderr, "                                      1: Use '<*>' notation to represent unobserved alleles\n");
    fprintf(stderr, "                                      2: Use '<NON_REF>' notation to represent unobserved alleles (a.k.a. GATK notation)\n");
    fprintf(stderr, "                                      3: Explode unobserved bases from {A,C,G,T} list\n");
    fprintf(stderr, "                                      4: Use '<*>' notation to represent unobserved alleles and explode unobserved bases from {A,C,G,T} list\n");
    fprintf(stderr, "                                      5: Use '<NON_REF>' notation to represent unobserved alleles and explode unobserved bases from {A,C,G,T} list\n");

    fprintf(stderr, "         -doGVCF [0]|1 ______________ 0: Disabled, 1: Output in gVCF format (requires: --rm-invar-sites 0, -doUnobserved 2, -addPL 1 and --gvcf-dps INT)\n");
    fprintf(stderr, "         -printPileup [0]|1 _________ 0: Disabled, 1: Also output in pileup format (<output_prefix>.pileup.gz)\n");
    fprintf(stderr, "         -printTruth [0]|1 __________ 0: Disabled, 1: Also output the VCF file containing the true genotypes (named <output_prefix>.truth.vcf)\n");
    fprintf(stderr, "         -printBasePickError [0]|1 __ 0: Disabled, 1: Print the base picking error probability to stdout.\n");
    fprintf(stderr, "                                      If --error-qs 1 is used, writes per-read base picking error probabilities to stdout.\n");
    fprintf(stderr, "                                      If --error-qs 0 or 2 is used, writes a single value which is used for all samples and sites.\n");
    fprintf(stderr, "         -printQsError [0]|1 ________ 0: Disabled, 1: Print the error probability used in quality score calculations to stdout.\n");
    fprintf(stderr, "                                      If --error-qs 2 is used, writes per-read quality score error probabilities to stdout.\n");
    fprintf(stderr, "                                      If --error-qs 0 or 1 is used, writes a single value which is used for all samples and sites.\n");
    fprintf(stderr, "         -printGlError [0]|1 ________ 0: Disabled, 1: Print the error probability used in genotype likelihood calculations to stdout. (requires: -GL 2)\n");
    fprintf(stderr, "                                      Since -GL 1 works directly with quality scores, this option is only available when -GL 2 is used.\n");
    fprintf(stderr, "                                      If --error-qs 2 is used, writes per-read error probabilities to stdout.\n");
    fprintf(stderr, "                                      If --error-qs 0 or 1 is used, writes a single value which is used for all samples and sites.\n");
    fprintf(stderr, "                                      If --precise-gl 1 is used, the printed values are the same as those printed by -printQsError.\n");
    fprintf(stderr, "         -printQScores [0]|1 ________ 0: Disabled, 1: Print the quality scores to stdout.\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "Output VCF/BCF tags:                  0: Do not add, 1: Add\n");
    fprintf(stderr, "         -addGL 0|[1] _______________ Genotype likelihoods (GL) tag\n");
    fprintf(stderr, "         -addGP [0]|1 _______________ Genotype probabilities (GP) tag\n");
    fprintf(stderr, "         -addPL [0]|1 _______________ Phred-scaled genotype likelihoods (PL) tag\n");
    fprintf(stderr, "         -addI16 [0]|1 ______________ I16 tag\n");
    fprintf(stderr, "         -addQS [0]|1 _______________ Quality score sum (QS) tag\n");
    fprintf(stderr, "         -addFormatDP [1]|0 _________ Per-sample read depth (FORMAT/DP) tag\n");
    fprintf(stderr, "         -addFormatAD [0]|1 _________ Per-sample allelic read depth (FORMAT/AD) tag\n");
    fprintf(stderr, "         -addFormatADF [0]|1 ________ Per-sample forward-strand allelic read depth (FORMAT/ADF) tag\n");
    fprintf(stderr, "         -addFormatADR [0]|1 ________ Per-sample reverse-strand allelic read depth (FORMAT/ADR) tag\n");
    fprintf(stderr, "         -addInfoDP [0]|1 ___________ Total read depth (INFO/DP) tag\n");
    fprintf(stderr, "         -addInfoAD [0]|1 ___________ Total allelic read depth (INFO/AD) tag\n");
    fprintf(stderr, "         -addInfoADF [0]|1 __________ Total forward-strand allelic read depth (INFO/ADF) tag\n");
    fprintf(stderr, "         -addInfoADR [0]|1 __________ Total reverse-strand allelic read depth (INFO/ADR) tag\n");


    fprintf(stderr, "\n");
}


/// @brief args_init - initialize the args struct and set arguments to default values
/// @return argStruct* - pointer to the args struct
argStruct* args_init() {

    args = (argStruct*)calloc(1, sizeof(argStruct));

    // -------------------------------------------- //
    // read in from command line

    args->verbose = 0;
    args->n_threads = 1;
    args->seed = -1;

    // --> Input/Output

    args->in_fn = NULL;
    args->gtSource = ARG_GTSOURCE_BINARY;
    args->out_fnprefix = NULL;
    args->output_mode = NULL;

    // --> Simulation parameters

    args->mps_depth = ARG_DEPTH_UNDEF;
    args->mps_depths_fn = NULL;
    args->error_rate = ARG_ERROR_RATE_UNDEF;
    args->error_qs = 0;
    args->beta_variance = -1.0;
    args->GL = 2;
    args->glModel1_theta = 0.83;
    args->platform = 0;
    args->usePreciseGlError = 0;
    args->i16_mapq = ARG_I16_MAPQ_DEFAULT;
    args->gvcf_dps_str = NULL;

    args->explode = 0;
    args->rmInvarSites = 0;
    args->rmEmptySites = 0;
    args->doUnobserved = ARG_DOUNOBSERVED_STAR;
    args->doGVCF = 0;
    args->printPileup = 0;
    args->printTruth = 0;
    args->printBasePickError = 0;
    args->printQsError = 0;
    args->printGlError = 0;
    args->printQScores = 0;

    args->addGL = 1;
    args->addGP = 0;
    args->addPL = 0;
    args->addI16 = 0;
    args->addQS = 0;
    args->addFormatDP = 1;
    args->addInfoDP = 0;
    args->addFormatAD = 0;
    args->addInfoAD = 0;
    args->addFormatADF = 0;
    args->addInfoADF = 0;
    args->addFormatADR = 0;
    args->addInfoADR = 0;


    // -------------------------------------------- //
    // set during program execution

    args->out_fn = NULL;
    args->out_truth_fn = NULL;
    args->out_pileup_fn = NULL;
    args->output_mode_str = NULL;

    args->in_fp = NULL;
    args->arg_fp = NULL;
    args->out_fp = NULL;
    args->out_truth_fp = NULL;
    args->out_pileup_fp = NULL;

    args->datetime = NULL;
    args->command = NULL;
    args->versionInfo = NULL;

    // -------------------------------------------- //
    // for internal use

    args->mps_depths = NULL;
    args->n_mps_depths = 1;
    args->betaSampler = NULL;
    args->poissonSampler = NULL;

    args->base_pick_error_prob = -1.0;

    args->preCalc = NULL;

    return args;
}


argStruct* args_get(int argc, char** argv) {

    args = args_init();

    if (0 == argc) {
        help_page();
        exit(0);
    }

    while (*argv) {
        char* arv = *argv;
        char* val = *(++argv);

        if ((strcmp("--help", arv) == 0) || (strcmp("-h", arv) == 0)) {
            help_page();
            exit(0);
        }

        else if ((strcmp("--version", arv) == 0) || (strcmp("-v", arv) == 0)) {
            version_page();
            exit(0);
        }

        else if ((strcmp("--verbose", arv) == 0) || (strcmp("-V", arv) == 0)) {
            args->verbose = atoi(val);
        }

        else if ((strcmp("--threads", arv) == 0) || (strcmp("-@", arv) == 0)) {
            args->n_threads = atoi(val);
        }

        else if ((strcmp("--seed", arv) == 0) || (strcmp("-s", arv) == 0)) {
            args->seed = atoi(val);
        }

        else if ((strcmp("--input", arv) == 0) || (strcmp("-i", arv) == 0)) {
            args->in_fn = strdup(val);
        }

        else if ((strcmp("--source", arv) == 0)) {
            args->gtSource = atoi(val);
        }

        else if ((strcmp("--output", arv) == 0) || (strcmp("-o", arv) == 0)) {
            args->out_fnprefix = strdup(val);
        }

        else if ((strcmp("--output-mode", arv) == 0) || (strcmp("-O", arv) == 0)) {
            args->output_mode = strdup(val);
        }

        else if ((strcmp("--depth", arv) == 0) || (strcmp("-d", arv) == 0)) {
            char* tmp = strdup(val);
            int isNumber = 1;
            int i = 0;
            if (tmp[0] == '-')
                i = 1;
            for (; tmp[i] != 0; ++i) {
                if (!isdigit(tmp[i])) {
                    if (tmp[i] == '.')
                        continue;
                    isNumber = 0;
                    if (strcasecmp("inf", tmp + i) == 0) {
                        args->mps_depth = ARG_DEPTH_INF;
                    }
                }
            }
            free(tmp);
            tmp = NULL;
            if (isNumber) {
                args->mps_depth = atof(val);
            } else if (ARG_DEPTH_INF != args->mps_depth) {
                NEVER;
            }
        }

        else if ((strcmp("--depths-file", arv) == 0) || (strcmp("-df", arv) == 0)) {
            args->mps_depths_fn = strdup(val);
            args->mps_depth = ARG_DEPTH_FILE;
        }

        else if ((strcmp("--error-rate", arv) == 0) || (strcmp("-e", arv) == 0)) {
            args->error_rate = atof(val);
        }

        else if ((strcasecmp("--error-qs", arv) == 0) || (strcasecmp("-eq", arv) == 0)) {
            args->error_qs = atoi(val);
        }

        else if ((strcasecmp("--beta-variance", arv) == 0) || (strcasecmp("-bv", arv) == 0)) {
            args->beta_variance = atof(val);
        }

        else if ((strcasecmp("--gl-model", arv) == 0) || (strcasecmp("-GL", arv) == 0)) {
            args->GL = atoi(val);
        }

        else if (strcasecmp("--gl1-theta", arv) == 0) {
            args->glModel1_theta = atof(val);
        }

        else if (strcasecmp("--platform", arv) == 0) {
            args->platform = atoi(val);
        }

        else if (strcasecmp("--precise-gl", arv) == 0) {
            args->usePreciseGlError = atoi(val);
        }

        else if (strcasecmp("--i16-mapq", arv) == 0) {
            args->i16_mapq = atoi(val);
        }

        else if (strcasecmp("--gvcf-dps", arv) == 0) {
            args->gvcf_dps_str = strdup(val);
        }

        else if (strcasecmp("-explode", arv) == 0) {
            args->explode = atoi(val);
        }

        else if (strcasecmp("--rm-invar-sites", arv) == 0) {
            args->rmInvarSites = atoi(val);
        }

        else if (strcasecmp("--rm-empty-sites", arv) == 0) {
            args->rmEmptySites = atoi(val);
        }

        else if (strcasecmp("-doUnobserved", arv) == 0) {
            args->doUnobserved = atoi(val);
        }

        else if (strcasecmp("-doGVCF", arv) == 0) {
            args->doGVCF = atoi(val);
        }

        else if (strcasecmp("-printPileup", arv) == 0) {
            args->printPileup = atoi(val);
        }

        else if (strcasecmp("-printTruth", arv) == 0) {
            args->printTruth = atoi(val);
        }

        else if (strcasecmp("-printBasePickError", arv) == 0) {
            args->printBasePickError = atoi(val);
        }

        else if (strcasecmp("-printQsError", arv) == 0) {
            args->printQsError = atoi(val);
        }

        else if (strcasecmp("-printGlError", arv) == 0) {
            args->printGlError = atoi(val);
        }

        else if (strcasecmp("-printQScores", arv) == 0) {
            args->printQScores = atoi(val);
        }

        else if ((strcasecmp("-addGL", arv) == 0) || (strcasecmp("-addFormatGL", arv) == 0)) {
            args->addGL = atoi(val);
        }

        else if ((strcasecmp("-addGP", arv) == 0) || (strcasecmp("-addFormatGP", arv) == 0)) {
            args->addGP = atoi(val);
        }

        else if ((strcasecmp("-addPL", arv) == 0) || (strcasecmp("-addFormatPL", arv) == 0)) {
            args->addPL = atoi(val);
        }

        else if ((strcasecmp("-addI16", arv) == 0) || (strcasecmp("-addFormatI16", arv) == 0)) {
            args->addI16 = atoi(val);
        }

        else if ((strcasecmp("-addQS", arv) == 0) || (strcasecmp("-addFormatQS", arv) == 0)) {
            args->addQS = atoi(val);
        }

        else if (strcasecmp("-addFormatDP", arv) == 0)
            args->addFormatDP = atoi(val);
        else if (strcasecmp("-addInfoDP", arv) == 0)
            args->addInfoDP = atoi(val);
        else if (strcasecmp("-addFormatAD", arv) == 0)
            args->addFormatAD = atoi(val);
        else if (strcasecmp("-addInfoAD", arv) == 0)
            args->addInfoAD = atoi(val);
        else if (strcasecmp("-addFormatADF", arv) == 0)
            args->addFormatADF = atoi(val);
        else if (strcasecmp("-addInfoADF", arv) == 0)
            args->addInfoADF = atoi(val);
        else if (strcasecmp("-addFormatADR", arv) == 0)
            args->addFormatADR = atoi(val);
        else if (strcasecmp("-addInfoADR", arv) == 0)
            args->addInfoADR = atoi(val);
        else {
            ERROR("Unknown arg:%s\n", arv);
        }
        ++argv;
    }

    // ---------------------------------------------------------------------- //
    // ARGUMENT VALIDATION

    CHECK_ARG_INTERVAL_INT(args->verbose, 0, 5, "--verbose");
    CHECK_ARG_INTERVAL_INT(args->n_threads, 0, ARG_NTHREADS_MAXVAL, "--threads");
    CHECK_ARG_INTERVAL_INT(args->seed, INT_MIN, INT_MAX, "--seed");

    // --> Input/Output

    // ** MUST BE SET **
    if (args->in_fn == NULL) {
        ERROR(
            "Input file is not specified. Please use -i/--input option to "
            "specify the input file.");
    }

    CHECK_ARG_INTERVAL_INT(args->gtSource, 0, 1, "--source");

    if (NULL == args->output_mode) {
        args->output_mode = strdup("b");
    }

    // --> Simulation parameters


    if (ARG_DEPTH_INF == args->mps_depth) {


        if (PROGRAM_WILL_SKIP_INPUT_HOMOREFGT_SITES) {
            ERROR("--rm-invar-sites %d cannot be used with -depth inf", args->rmInvarSites);
        }
        if (PROGRAM_WILL_SKIP_INPUT_HOMOALTGT_SITES) {
            ERROR("--rm-invar-sites %d cannot be used with -depth inf", args->rmInvarSites);
        }
        if (args->printPileup != 0) {
            ERROR(
                "(-printPileup 1) Cannot output pileup file when --depth inf is set.");
        }
        if (args->addFormatDP) {
            ERROR(
                "(-addFormatDP 1) FORMAT/DP tag cannot be added when --depth "
                "inf is set.");
        }

        if (args->addInfoDP) {
            ERROR(
                "(-addInfoDP 1) INFO/DP tag cannot be added when --depth inf "
                "is set.");
        }

        if (args->addI16 != 0) {
            ERROR(
                "(-addI16 != 0) I16 tag cannot be added when --depth inf is set.");
        }

        if (1 == args->addQS) {
            ERROR("(-addQS 1) QS tag cannot be added when --depth inf is set.");
        }

        if (args->addFormatAD) {
            ERROR(
                "(-addFormatAD 1) FORMAT/AD tag cannot be added when --depth "
                "inf is set.");
        }
        if (args->addFormatADF) {
            ERROR(
                "(-addFormatADF 1) FORMAT/ADF tag cannot be added when --depth "
                "inf is set.");
        }
        if (args->addFormatADR) {
            ERROR(
                "(-addFormatADR 1) FORMAT/ADR tag cannot be added when --depth "
                "inf is set.");
        }

        if (args->addInfoAD) {
            ERROR(
                "(-addInfoAD 1) INFO/AD tag cannot be added when --depth inf "
                "is set.");
        }
        if (args->addInfoADF) {
            ERROR(
                "(-addInfoADF 1) INFO/ADF tag cannot be added when --depth inf "
                "is set.");
        }
        if (args->addInfoADR) {
            ERROR(
                "(-addInfoADR 1) INFO/ADR tag cannot be added when --depth inf "
                "is set.");
        }

        if (0 != args->error_rate) {
            ERROR(
                "Cannot simulate true values (--depth inf) with error rate "
                "(--error-rate) set to %f. Please set --error-rate to 0 and "
                "rerun.",
                args->error_rate);
        }
    } else if (ARG_DEPTH_UNDEF == args->mps_depth) {
        ERROR("Average per-site read depth value is required. Please set it using --depth or --depths-file and re-run.");
    } else if (ARG_DEPTH_FILE == args->mps_depth) {
        //
    } else {
        CHECK_ARG_INTERVAL_DBL(args->mps_depth, 0.0, ARG_DEPTH_MAXVAL, "--depth");
    }


    if (args->error_rate == ARG_ERROR_RATE_UNDEF) {
        ERROR("Error rate is not specified. Please use --error-rate option to specify the error rate. Allowed range: [0.0, 1.0]");
    }
    CHECK_ARG_INTERVAL_IE_DBL(args->error_rate, 0.0, 1.0, "--error-rate");

    CHECK_ARG_INTERVAL_INT(args->error_qs, 0, 2, "--error-qs");

    if (!PROGRAM_BETA_VAR_IS_UNDEF) {
        CHECK_ARG_INTERVAL_DBL(args->beta_variance, 0.0, 1.0, "--beta-variance");
    }

    CHECK_ARG_INTERVAL_INT(args->GL, 1, 2, "--gl-model");
    CHECK_ARG_INTERVAL_DBL(args->glModel1_theta, 0.0, 1.0, "--gl1-theta");
    CHECK_ARG_INTERVAL_01(args->platform, "--platform");
    CHECK_ARG_INTERVAL_01(args->usePreciseGlError, "--precise-gl");
    CHECK_ARG_INTERVAL_INT(args->i16_mapq, 0, 60, "--i16-mapq");

    CHECK_ARG_INTERVAL_01(args->explode, "-explode");
    CHECK_ARG_INTERVAL_INT(args->rmInvarSites, 0, 7, "--rm-invar-sites");
    CHECK_ARG_INTERVAL_01(args->rmEmptySites, "--rm-empty-sites");
    CHECK_ARG_INTERVAL_INT(args->doUnobserved, 0, 5, "-doUnobserved");
    CHECK_ARG_INTERVAL_01(args->doGVCF, "-doGVCF");
    CHECK_ARG_INTERVAL_01(args->printPileup, "-printPileup");
    CHECK_ARG_INTERVAL_01(args->printTruth, "-printTruth");
    CHECK_ARG_INTERVAL_01(args->printBasePickError, "-printBasePickError");
    CHECK_ARG_INTERVAL_01(args->printQsError, "-printQsError");
    CHECK_ARG_INTERVAL_01(args->printGlError, "-printGlError");
    CHECK_ARG_INTERVAL_01(args->printQScores, "-printQScores");

    CHECK_ARG_INTERVAL_01(args->addGL, "-addGL");
    CHECK_ARG_INTERVAL_01(args->addGP, "-addGP");
    CHECK_ARG_INTERVAL_01(args->addPL, "-addPL");
    CHECK_ARG_INTERVAL_INT(args->addI16, 0, 2, "-addI16");
    CHECK_ARG_INTERVAL_01(args->addQS, "-addQS");
    CHECK_ARG_INTERVAL_01(args->addFormatDP, "-addFormatDP");
    CHECK_ARG_INTERVAL_01(args->addInfoDP, "-addInfoDP");
    CHECK_ARG_INTERVAL_01(args->addFormatAD, "-addFormatAD");
    CHECK_ARG_INTERVAL_01(args->addInfoAD, "-addInfoAD");
    CHECK_ARG_INTERVAL_01(args->addFormatADF, "-addFormatADF");
    CHECK_ARG_INTERVAL_01(args->addInfoADF, "-addInfoADF");
    CHECK_ARG_INTERVAL_01(args->addFormatADR, "-addFormatADR");
    CHECK_ARG_INTERVAL_01(args->addInfoADR, "-addInfoADR");

    // ---------------------------------------------------------------------- //
    // ARGUMENT COMBINATION COMPATIBILITY CHECKS

    if (args->usePreciseGlError && (1 == args->GL)) {
        ERROR("Precise genotype likelihood error (--precise-gl 1) is not supported with genotype likelihood model 1 (--gl-model 1).");
    }

    if ((args->beta_variance >= 0) && (0 == args->error_qs)) {
        ERROR("--beta-variance %e requires --error-qs 1 or 2.", args->beta_variance);
    }

    if (1 == args->error_qs || 2 == args->error_qs) {
        if (0.0 == args->error_rate) {
            ERROR("--error-qs 1 or 2 requires --error-rate > 0 (found %f). If you want to simulate reads without errors, please use --error-qs 0 instead.", args->error_rate);
        } else if (args->error_rate > 0.0) {
            if (0.0 == args->beta_variance) {
                ERROR("--error-qs 1 or 2 requires --beta-variance > 0 (found %e). If you want to simulate quality scores without errors, please use --error-qs 0 instead.", args->beta_variance);
            } else if (args->beta_variance > 0.0) {
            } else {
                ERROR("--beta-variance should be set to a positive value. Please set --beta-variance and rerun.");
            }

        } else {
            ERROR("Error rate (--error-rate) should be set to a positive value. Please set --error-rate and rerun.");
        }
    }


    if (1 == args->doGVCF) {

        if (!args->addFormatDP) {
            ERROR("[-doGVCF 1] -addFormatDP 1 is required for gVCF output. Please set -addFormatDP 1 and rerun.");
        }

        WARN("\n-> [-doGVCF 1] GVCF mode is under development. Please use with caution and report any issues.");

        if (0 != args->rmInvarSites) {
            ERROR("\n-> [-doGVCF 1] --rm-invar-sites 0 is required. Please set --rm-invar-sites 0 and rerun.");
        }

        if (NULL == args->gvcf_dps_str) {
            ERROR("\n-> [-doGVCF 1] --gvcf-dps is required. Please set --gvcf-dps and rerun.");
        }

        if (!PROGRAM_WILL_ADD_UNOBSERVED) {
            ERROR("\n-> [-doGVCF 1] Adding unobserved alleles is required for gVCF output. Please set -doUnobserved to %d or %d and rerun.", ARG_DOUNOBSERVED_STAR, ARG_DOUNOBSERVED_NONREF);
        }

        if (ARG_DEPTH_INF == args->mps_depth) {
            ERROR("\n-> [-doGVCF 1] --depth inf is not supported with gVCF output. Please set --depth to a finite value and rerun.");
        }

    } else {
        if (args->gvcf_dps_str != NULL) {
            ERROR("\n-> [--gvcf-dps] --gvcf-dps requires -doGVCF 1. Please set -doGVCF 1 and rerun.");
        }
    }


    if (args->printGlError && (1 == args->GL)) {
        ERROR("\n-> [-printGlError 1] Printing the error probability used in genotype likelihood calculations (-printGlError 1) is not supported with genotype likelihood model 1 (--gl-model 1).");
    }


    if ((0 == args->addI16) && (ARG_I16_MAPQ_DEFAULT != args->i16_mapq)) {
        ERROR("\n-> [--i16-mapq] --i16-mapq is set, but printing I16 tag is disabled. Please set -addI16 to 1, or remove --i16-mapq and rerun.");
    }

    // ---------------------------------------------------------------------- //
    // SET ARGUMENT VARIABLES

    if (NULL == args->out_fnprefix) {
        WARN("Output file prefix is not specified. Setting it to 'output'.\n");
        args->out_fnprefix = strdup("output");
    }

    if (PROGRAM_WILL_ADD_STAR) {
        nonref_str = "<*>";
    } else if (PROGRAM_WILL_ADD_NONREF) {
        nonref_str = "<NON_REF>";
    }


    if (ARG_DEPTH_INF == args->mps_depth) {
        if (PROGRAM_WILL_SKIP_SIM_INVAR_SITES) {
            ERROR("[--rm-invar-sites %d] Cannot skip invariable sites when --depth inf is set.", args->rmInvarSites);
        }
        if (1 == args->doGVCF) {
            ERROR("[-doGVCF 1] Cannot output gVCF when --depth inf is set.");
        }
    }

    if (ARG_DEPTH_FILE == args->mps_depth) {
        args->mps_depths = args->read_depthsFile();
    }

    args->arg_fp = open_FILE(args->out_fnprefix, ".arg");

    args->datetime = strdup(get_time());

    if (args->error_qs != 0) {

#if __USE_STD_BETA__ == 1
        args->betaSampler = new BetaSampler(args->error_rate, args->beta_variance, args->seed);
#else
        args->betaSampler = BetaSampler_init(args->error_rate, args->beta_variance, args->seed);
#endif
    }

    // -> SEED

    if (args->seed == -1) {
        int rseed = time(NULL);
        fprintf(stderr, "\n-> No seed was given. Setting the random seed to the randomly chosen value: %d\n", rseed);
        args->seed = rseed;
    }

    // set seed for rng0
    srand48(args->seed);

    // set seed for rng1 
    rng1_seeder[1] = (unsigned short)((long)args->seed);
    rng1_seeder[2] = (unsigned short)(((long)args->seed) >> 16);

    rng2_seeder[1] = (unsigned short)((long)args->seed);
    rng2_seeder[2] = (unsigned short)(((long)args->seed) >> 16);

    if (ARG_DEPTH_INF != args->mps_depth) {
        if (1 == args->n_mps_depths) {
            args->poissonSampler = (PoissonSampler**)malloc(sizeof(PoissonSampler*));
            args->poissonSampler[0] = PoissonSampler_init(args->mps_depth);
        } else {
            //array of poisson samplers
            args->poissonSampler = (PoissonSampler**)malloc(args->n_mps_depths * sizeof(PoissonSampler*));
            for (int i = 0; i < args->n_mps_depths; ++i) {
                args->poissonSampler[i] = PoissonSampler_init(args->mps_depths[i]);
            }
        }
    }


    if (args->printBasePickError) {
        fprintf(stderr, "[-printBasePickError 1] Program will print the base picking error probability to stdout.\n");
    }

    if (args->printGlError) {
        fprintf(stderr, "[-printGlError 1] Program will print the error probability used in genotype likelihood calculations to stdout.\n");
    }

    if (args->printQScores) {
        fprintf(stderr, "[-printQScores 1] Program will print the quality scores to stdout.\n");
    }

    if (0 == args->error_qs) {
        args->base_pick_error_prob = args->error_rate;
        if (args->printBasePickError) {
            fprintf(stdout, "base_pick_error_prob\tNA\tNA\tNA\tNA\t%f\n", args->base_pick_error_prob);
        }
    } else if (1 == args->error_qs) {
        args->base_pick_error_prob = -1.0; // sampled
    } else if (2 == args->error_qs) {
        args->base_pick_error_prob = args->error_rate;
        if (args->printBasePickError) {
            fprintf(stdout, "base_pick_error_prob\tNA\tNA\tNA\tNA\t%f\n", args->base_pick_error_prob);
        }
    }


    // if ((0 == args->error_qs) || (1 == args->error_qs)) {
        // args->preCalc = new preCalcStruct();
        // args->preCalc->error_prob_forQs = args->error_rate;
        // args->preCalc->qScore = get_qScore(args->preCalc->error_prob_forQs);
        // args->preCalc->q5 = args->preCalc->qScore << 5;
        // if (1 == args->GL) {
            // glModel1 = glModel1_init();
            // calculate_gls = alleles_calculate_gls_log10_glModel1_fixedQScore;
        // } else if (2 == args->GL) {
            // args->preCalc->prepare_gls_preCalc();
            // calculate_gls = alleles_calculate_gls_log10_glModel2_fixedQScore;
        // }
//
        // if (args->printGlError) {
            // fprintf(stdout, "gl_error_prob\tNA\tNA\tNA\tNA\t%f\n", args->preCalc->error_prob_forGl);
        // }
//
        // if (args->printQsError) {
            // fprintf(stdout, "qs_error_prob\tNA\tNA\tNA\tNA\t%f\n", args->preCalc->error_prob_forQs);
        // }
//
        // if (args->printQScores) {
            // fprintf(stdout, "qs\tNA\tNA\tNA\tNA\t%d\n", args->preCalc->qScore);
        // }
//
    // } else if (2 == args->error_qs) {
        // if (1 == args->GL) {
            // glModel1 = glModel1_init();
            // calculate_gls = alleles_calculate_gls_log10_glModel1;
        // } else if (2 == args->GL) {
            // if (args->usePreciseGlError) {
                // calculate_gls = alleles_calculate_gls_log10_glModel2_precise1;
            // } else {
                // calculate_gls = alleles_calculate_gls_log10_glModel2_precise0;
            // }
        // }
    // }
//
//
    // ---------------------------------------------------------------------- //
    // PRINT ARGUMENTS

    ASSERT(asprintf(&args->versionInfo, "vcfgl [version: %s] [build: %s %s] [htslib: %s]\n", VCFGL_VERSION, __DATE__, __TIME__, hts_version()));
    fprintf(stderr, "\n%s\n", args->versionInfo);
    fprintf(args->arg_fp, args->versionInfo);

    char depth_val[1024];
    if (ARG_DEPTH_INF == args->mps_depth) {
        sprintf(depth_val, "--depth %s", "inf");
    } else if (ARG_DEPTH_FILE == args->mps_depth) {
        ASSERT(NULL != args->mps_depths_fn);
        sprintf(depth_val, "--depths-file %s", args->mps_depths_fn);
    } else {
        sprintf(depth_val, "--depth %f", args->mps_depth);
    }

    char gvcf_dps_val[1024];
    if (NULL == args->gvcf_dps_str) {
        gvcf_dps_val[0] = '\0';
    } else {
        sprintf(gvcf_dps_val, "--gvcf-dps %s", args->gvcf_dps_str);
    }

    char gl_model_str[256];
    if (1 == args->GL) {
        sprintf(gl_model_str, "--gl-model %d --gl1-theta %f", args->GL, args->glModel1_theta);
    } else {
        sprintf(gl_model_str, "--gl-model %d", args->GL);
    }

    char beta_variance_str[256];
    if (args->beta_variance >= 0) {
        sprintf(beta_variance_str, "--beta-variance %e", args->beta_variance);
    } else {
        beta_variance_str[0] = '\0';
    }

    char i16_mapq_str[256];
    if (args->addI16 == 1) {
        sprintf(i16_mapq_str, "--i16-mapq %d", args->i16_mapq);
    } else {
        i16_mapq_str[0] = '\0';
    }

    ASSERT(asprintf(
        &args->command,
        "Command: vcfgl --verbose %d --threads %d --seed %d --input %s --source %d --output %s --output-mode %s %s --error-rate %f --error-qs %d %s %s --platform %d --precise-gl %d --i16-mapq %d %s -explode %d --rm-invar-sites %d --rm-empty-sites %d -doUnobserved %d -doGVCF %d -printPileup %d -printTruth %d  -printBasePickError %d -printQsError %d -printGlError %d -printQScores %d -addGL %d -addGP %d -addPL %d -addI16 %d -addQS %d -addFormatDP %d -addInfoDP %d -addFormatAD %d -addInfoAD %d -addFormatADF %d -addInfoADF %d -addFormatADR %d -addInfoADR %d",
        args->verbose,
        args->n_threads,
        args->seed,
        args->in_fn,
        args->gtSource,
        args->out_fnprefix,
        args->output_mode,
        depth_val,
        args->error_rate,
        args->error_qs,
        beta_variance_str,
        gl_model_str,
        args->platform,
        args->usePreciseGlError,
        args->i16_mapq,
        gvcf_dps_val,
        args->explode,
        args->rmInvarSites,
        args->rmEmptySites,
        args->doUnobserved,
        args->doGVCF,
        args->printPileup,
        args->printTruth,
        args->printBasePickError,
        args->printQsError,
        args->printGlError,
        args->printQScores,
        args->addGL,
        args->addGP,
        args->addPL,
        args->addI16,
        args->addQS,
        args->addFormatDP,
        args->addInfoDP,
        args->addFormatAD,
        args->addInfoAD,
        args->addFormatADF,
        args->addInfoADF,
        args->addFormatADR,
        args->addInfoADR) > 0);



    fprintf(stderr, "\n%s\n\n", args->command);
    fprintf(stderr, "\n%s\n", args->datetime);
    fprintf(args->arg_fp, "\n%s\n\n", args->command);
    fprintf(args->arg_fp, "\n%s\n", args->datetime);

    if (args->printPileup) {
        args->out_pileup_fn = (char*)malloc(strlen(args->out_fnprefix) + strlen(".pileup.gz") + 1);
        strcpy(args->out_pileup_fn, args->out_fnprefix);
        strcat(args->out_pileup_fn, ".pileup.gz");
    }


    switch (*args->output_mode) {
    case 'v':
        if (args->n_threads > 1) {
            ERROR(
                "Multithreading is not supported for VCF output. Please "
                "set --threads 1 and rerun.");
        }
        args->out_fn = (char*)malloc(strlen(args->out_fnprefix) + strlen(".vcf") + 1);
        strcpy(args->out_fn, args->out_fnprefix);
        strcat(args->out_fn, ".vcf");
        args->output_mode_str = strdup("w");

        if (args->printTruth) {
            args->out_truth_fn = (char*)malloc(strlen(args->out_fnprefix) + strlen(".truth.vcf") + 1);
            strcpy(args->out_truth_fn, args->out_fnprefix);
            strcat(args->out_truth_fn, ".truth.vcf");
        }
        break;
    case 'b':
        args->out_fn =
            (char*)malloc(strlen(args->out_fnprefix) + strlen(".bcf") + 1);
        strcpy(args->out_fn, args->out_fnprefix);
        strcat(args->out_fn, ".bcf");
        args->output_mode_str = strdup("wb");

        if (args->printTruth) {
            args->out_truth_fn = (char*)malloc(strlen(args->out_fnprefix) + strlen(".truth.bcf") + 1);
            strcpy(args->out_truth_fn, args->out_fnprefix);
            strcat(args->out_truth_fn, ".truth.bcf");
        }
        break;
    case 'z':
        if (args->n_threads > 1) {
            ERROR(
                "Multithreading is not supported for VCF output. Please "
                "set --threads 1 and rerun.");
        }
        args->out_fn =
            (char*)malloc(strlen(args->out_fnprefix) + strlen(".vcf.gz") + 1);
        strcpy(args->out_fn, args->out_fnprefix);
        strcat(args->out_fn, ".vcf.gz");
        args->output_mode_str = strdup("wz");

        if (args->printTruth) {
            args->out_truth_fn = (char*)malloc(strlen(args->out_fnprefix) + strlen(".truth.vcf.gz") + 1);
            strcpy(args->out_truth_fn, args->out_fnprefix);
            strcat(args->out_truth_fn, ".truth.vcf.gz");
        }
        break;
    case 'u':
        args->out_fn =
            (char*)malloc(strlen(args->out_fnprefix) + strlen(".bcf") + 1);
        strcpy(args->out_fn, args->out_fnprefix);
        strcat(args->out_fn, ".bcf");
        args->output_mode_str = strdup("wbu");

        if (args->printTruth) {
            args->out_truth_fn = (char*)malloc(strlen(args->out_fnprefix) + strlen(".truth.bcf") + 1);
            strcpy(args->out_truth_fn, args->out_fnprefix);
            strcat(args->out_truth_fn, ".truth.bcf");
        }
        break;
    }

    fprintf(stderr, "\n-> Log file: %s.arg\n", args->out_fnprefix);
    fprintf(args->arg_fp, "\n-> Log file: %s.arg\n", args->out_fnprefix);
    fprintf(stderr, "-> Simulation output file: %s\n", args->out_fn);
    fprintf(args->arg_fp, "-> Simulation output file: %s\n", args->out_fn);
    if (args->printTruth) {
        fprintf(stderr, "-> True genotypes output file: %s\n", args->out_truth_fn);
        fprintf(args->arg_fp, "-> True genotypes output file: %s\n", args->out_truth_fn);
    }
    if (args->printPileup) {
        fprintf(stderr, "-> Pileup output file: %s\n", args->out_pileup_fn);
        fprintf(args->arg_fp, "-> Pileup output file: %s\n", args->out_pileup_fn);
    }


    // ---------------------------------------------------------------------- //
    // ARGUMENT WARNINGS

    if (args->doUnobserved == ARG_DOUNOBSERVED_TRIM) {
        WARN("[-doUnobserved %d] Unobserved alleles will be trimmed from the output. Please note that this is for experimental purposes only, and the invariable sites may not be accepted by some downstream tools due to the lack of information about the unobserved alleles.", args->doUnobserved);
    }

    if (1 == args->addI16) {
        WARN("[-addI16 1] This feature is experimental. Please use with caution and report any issues.");
    }

    if (1 == args->error_qs) {
        // Per-site error probabilities will be sampled from a beta distribution. 
        WARN("[--error-qs 1] This feature is experimental. Please use with caution and report any issues.");
    }

    if (1 != args->n_threads) {
        if (0 == args->n_threads) {
            VWARN("--threads 0 implies --threads 1. Setting number of threads to 1 (no multithreading).");
            args->n_threads = 1;
        }
    }

    if (ARG_DEPTH_INF == args->mps_depth) {
        WARN("[--depth inf] This feature is experimental. Please use with caution and report any issues.");
    }

    return (args);
}

void args_destroy(argStruct* args) {

    // --> Input/Output

    free(args->in_fn);
    args->in_fn = NULL;
    free(args->out_fnprefix);
    args->out_fnprefix = NULL;
    free(args->output_mode);
    args->output_mode = NULL;

    if (NULL != args->mps_depths_fn) {
        free(args->mps_depths_fn);
        args->mps_depths_fn = NULL;
    }

    free(args->gvcf_dps_str);
    args->gvcf_dps_str = NULL;

    // -------------------------------------------- //
    // set during program execution

    free(args->out_fn);
    args->out_fn = NULL;
    if (NULL != args->out_truth_fn) {
        free(args->out_truth_fn);
        args->out_truth_fn = NULL;
    }
    if (NULL != args->out_pileup_fn) {
        free(args->out_pileup_fn);
        args->out_pileup_fn = NULL;
    }
    free(args->output_mode_str);
    args->output_mode_str = NULL;

    ASSERT(0 == fclose(args->arg_fp));
    args->arg_fp = NULL;

    free(args->datetime);
    args->datetime = NULL;
    free(args->command);
    args->command = NULL;
    free(args->versionInfo);
    args->versionInfo = NULL;

    // -------------------------------------------- //
    // for internal use
    if (NULL != args->mps_depths) {
        free(args->mps_depths);
        args->mps_depths = NULL;
    }

    if (args->betaSampler != NULL) {

#if __USE_STD_BETA__==1
        delete (args->betaSampler);
#else
        BetaSampler_destroy(args->betaSampler);
#endif
    }
    if (args->poissonSampler != NULL) {
        for (int i = 0;i < args->n_mps_depths;++i) {
            free(args->poissonSampler[i]);
            args->poissonSampler[i] = NULL;
        }
        free(args->poissonSampler);
        args->poissonSampler = NULL;
    }

    if (args->preCalc != NULL) {
        delete (args->preCalc);
    }

    free(args);
    args = NULL;
}

size_t fsize(const char* fname) {
    struct stat st;
    stat(fname, &st);
    return st.st_size;
}


/// @brief read depths file args->mps_depths_fn into args->mps_depths array and set args->n_mps_depths with array size
double* argStruct::read_depthsFile(void) {

    if (this->verbose > 0) {
        fprintf(stderr, "-> Reading depths file: %s\n", this->mps_depths_fn);
    }

    FILE* fp = NULL;
    if ((fp = fopen(this->mps_depths_fn, "r")) == NULL) {
        ERROR("Could not open file: %s\n", this->mps_depths_fn);
    }

    double val;
    char* lineend;

    int n = 0;
    char buf[255] = { '\0' };
    int ntmp = BUFSIZE_NINDS;
    double* ret = (double*)malloc(ntmp * sizeof(double));

    while (fgets(buf, 255, fp)) {

        val = strtod(buf, &lineend);
        ASSERT(lineend != buf);

        if (n == ntmp) {
            ntmp *= 2;
            ret = (double*)realloc(ret, ntmp * sizeof(double));
        }

        ret[n] = val;
        ++n;

    }

    if (n == 0) {
        ERROR("Could not read any values from depths file: %s\n", this->mps_depths_fn);
    }

    if (n != ntmp) {
        ret = (double*)realloc(ret, n * sizeof(double));
    }
    this->n_mps_depths = n;

    if (this->verbose > 0) {
        fprintf(stderr, "\n-> Read %d values from depths file: %s\n", n, this->mps_depths_fn);

        if (this->verbose > 1) {
            fprintf(stderr, "\t-> Values: ");
            for (int i = 0; i < n; ++i) {
                fprintf(stderr, "%f ", ret[i]);
            }
            fprintf(stderr, "\n");
        }
    }

    fclose(fp);
    return(ret);
}


