#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/vcf.h>
#include <math.h>
#include <stdbool.h>

#define SETALLELES_VERSION "0.0.1"

// === ANSI codes ===
#define ANSI_ESC          "\033"
#define ANSI_RESET        ANSI_ESC "[0m"
#define ANSI_BOLD_FMT     ANSI_ESC "[1m"
#define ANSI_COLOR_FMT    ANSI_ESC "[1;%dm"

// === Color values ===
#define ANSI_COLOR_BLACK    30
#define ANSI_COLOR_RED      31
#define ANSI_COLOR_GREEN    32
#define ANSI_COLOR_YELLOW   33
#define ANSI_COLOR_BLUE     34

// === Color name aliases (what you want to use) ===
#define BLACK    ANSI_COLOR_BLACK
#define RED      ANSI_COLOR_RED
#define GREEN    ANSI_COLOR_GREEN
#define YELLOW   ANSI_COLOR_YELLOW
#define BLUE     ANSI_COLOR_BLUE

/*
* Macro:[FPRINTF_BOLD_COLOR]
* print a custom message in bold color
*/
#define FPRINTF_BOLD_COLOR(stream, color, ...)                      \
    do {                                                            \
        fprintf((stream), ANSI_COLOR_FMT, (color));                 \
        fprintf((stream), __VA_ARGS__);                             \
        fprintf((stream), ANSI_RESET);                              \
    } while (0)

/*
* Macro:[FPRINTF_BOLD]
* print a custom message in bold
*/
#define FPRINTF_BOLD(stream, ...)                                       \
    do {                                                                \
        fprintf((stream), ANSI_BOLD_FMT);                               \
        fprintf((stream), __VA_ARGS__);                                 \
    } while (0)

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
        fprintf(stderr, "\n\n*******\n[ERROR](%s)<%s:%d>\n\t", __func__, \
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
                    __FILE__, __func__, __LINE__, #expr);               \
            exit(1);                                                        \
        }                                                                   \
    } while (0);

/*
* Macro:[WARN]
* print a custom warning message
*/
#define WARN(...)                                                            \
    do {                                                                     \
        fprintf(stderr, "\n\n[WARNING](%s/%s:%d): ", __FILE__, __func__, \
                __LINE__);                                                   \
        fprintf(stderr, __VA_ARGS__);                                        \
        fprintf(stderr, "\n");                                               \
    } while (0);


#define MAX_LINE 256
#define MAX_NALTS 4


typedef struct {
    // ref[n_sites] = char* ref
    char** ref;
    // alts[n_sites][n_alts] = char* alt
    char*** alts;
    // n_alts[n_sites] = uint8_t n_alts
    uint8_t* n_alts;
    uint64_t n_sites;
} alleles_t;

alleles_t *alleles_init(void){
    alleles_t *alleles = NULL;
    alleles = (alleles_t*) malloc(sizeof(alleles_t));
    ASSERT(alleles!=NULL);
    alleles->ref = NULL;
    alleles->alts = NULL;
    alleles->n_alts = NULL;
    alleles->n_sites = 0;
    return(alleles);
}

void alleles_destroy(alleles_t *alleles){
    if (alleles == NULL) return;
    
    for(uint64_t i = 0; i < alleles->n_sites; ++i){
        if (alleles->ref != NULL && alleles->ref[i] != NULL) {
            free(alleles->ref[i]);
            alleles->ref[i] = NULL;
        }
        if (alleles->alts != NULL && alleles->alts[i] != NULL) {
            for(uint8_t j = 0; j < alleles->n_alts[i]; ++j){
                if (alleles->alts[i][j] != NULL) {
                    free(alleles->alts[i][j]);
                    alleles->alts[i][j] = NULL;
                }
            }
            free(alleles->alts[i]);
            alleles->alts[i] = NULL;
        }
    }
    if (alleles->ref != NULL) {
        free(alleles->ref);
        alleles->ref = NULL;
    }
    if (alleles->alts != NULL) {
        free(alleles->alts);
        alleles->alts = NULL;
    }
    if (alleles->n_alts != NULL) {
        free(alleles->n_alts);
        alleles->n_alts = NULL;
    }
    free(alleles);
    alleles = NULL;
}

void alleles_add_site(alleles_t *alleles, char* ref, char** alts, uint8_t n_alts){
    alleles->ref = (char**) realloc(alleles->ref, (alleles->n_sites + 1) * sizeof(char*));
    ASSERT(alleles->ref!=NULL);
    alleles->n_alts = (uint8_t*) realloc(alleles->n_alts, (alleles->n_sites + 1) * sizeof(uint8_t));
    ASSERT(alleles->n_alts!=NULL);
    alleles->n_alts[alleles->n_sites] = n_alts;
    alleles->alts = (char***) realloc(alleles->alts, (alleles->n_sites + 1) * sizeof(char**));
    ASSERT(alleles->alts!=NULL);
    alleles->alts[alleles->n_sites] = (char**) malloc(n_alts * sizeof(char*));
    ASSERT(alleles->alts[alleles->n_sites]!=NULL);

    alleles->ref[alleles->n_sites] = strdup(ref);   
    ASSERT(alleles->ref[alleles->n_sites]!=NULL);
    for(uint8_t i = 0; i < n_alts; ++i){
        alleles->alts[alleles->n_sites][i] = strdup(alts[i]);
        ASSERT(alleles->alts[alleles->n_sites][i]!=NULL);
    }
    alleles->n_sites++;
}

alleles_t *load_alleles_tsv(const char *fname) {
    FILE *fp = fopen(fname, "r");
    ASSERT(fp!=NULL);

    alleles_t *alleles = alleles_init();
    char line[MAX_LINE];

    char* ref = NULL;
    char* altsline = NULL;
    char* alts[MAX_NALTS] = {NULL};
    uint8_t n_alts = 0;

    while (fgets(line, MAX_LINE, fp)) {

        ref = strtok(line, "\t\n");
        ASSERT(ref!=NULL);
        altsline = strtok(NULL, "\t\n");
        ASSERT(altsline!=NULL);

        char *alt_tok = strtok(altsline, ",");
        ASSERT(alt_tok!=NULL);
        while(alt_tok != NULL){
            alts[n_alts++] = strdup(alt_tok);
            alt_tok = strtok(NULL, ",");
        }

        alleles_add_site(alleles, ref, alts, n_alts);

        // clean up
        ref = NULL;
        altsline = NULL;
        for(uint8_t i = 0; i < n_alts; ++i){
            free(alts[i]);
            alts[i] = NULL;
        }
        n_alts = 0;
    }

    //// print alleles
    //printf("\n\n");
    //for(uint64_t i = 0; i < alleles->n_sites; ++i){
    //    printf("ref: %s\n", alleles->ref[i]);
    //    for(uint8_t j = 0; j < alleles->n_alts[i]; ++j){
    //        printf("alts[%u]: %s\n", j, alleles->alts[i][j]);
    //    }
    //}
    //NEVER;

    ASSERT(fclose(fp) == 0);
    return alleles;
}

typedef struct{
    char* input_bcf_file;
    char* alleles_tsv_file;
    char* output_mode;
    char* output_prefix;
    char* output_fn;
} args_t;

args_t* args_init(void){
    args_t* args = (args_t*) malloc(sizeof(args_t));
    ASSERT(args!=NULL);
    args->input_bcf_file = NULL;
    args->alleles_tsv_file = NULL;
    args->output_mode = NULL;
    args->output_prefix = NULL;
    args->output_fn = NULL;
    return(args);
}

void args_destroy(args_t* args){
    if (args == NULL) return;
    
    if (args->input_bcf_file != NULL) {
        free(args->input_bcf_file);
        args->input_bcf_file = NULL;
    }
    if (args->alleles_tsv_file != NULL) {
        free(args->alleles_tsv_file);
        args->alleles_tsv_file = NULL;
    }
    if (args->output_mode != NULL) {
        free(args->output_mode);
        args->output_mode = NULL;
    }
    if (args->output_prefix != NULL) {
        free(args->output_prefix);
        args->output_prefix = NULL;
    }
    if (args->output_fn != NULL) {
        free(args->output_fn);
        args->output_fn = NULL;
    }
    free(args);
    args = NULL;
}

void set_arg_value(char** arg_to_set, const char* arg_str, const char* val) {
    ASSERT(arg_to_set!=NULL);
    if (NULL == val || strlen(val) == 0) {
        ERROR("Argument '%s' requires a value.", arg_str);
    }
    if (*arg_to_set != NULL) {
        free(*arg_to_set);
        *arg_to_set = NULL;
    }
    *arg_to_set = strdup(val);
    ASSERT(*arg_to_set!=NULL);
}
    

void version_page(void) {
    fprintf(stderr, "setAlleles [version: %s] [build: %s %s] [htslib: %s]\n", SETALLELES_VERSION, __DATE__, __TIME__, hts_version());
    fprintf(stderr, "\n");
}

void help_page(void){

    fprintf(stderr, "\n");
    FPRINTF_BOLD(stderr, "----------------------------------------------------------\n");
    FPRINTF_BOLD(stderr, "Program: setAlleles\n");
    FPRINTF_BOLD(stderr, "License: GNU GPLv3.0\n");
    FPRINTF_BOLD(stderr, "Version: %s (htslib: %s)\n", SETALLELES_VERSION, hts_version());
    FPRINTF_BOLD(stderr, "Build: %s %s\n", __DATE__, __TIME__);
    FPRINTF_BOLD(stderr, "-> setAlleles is a subprogram of vcfgl.\n");
    FPRINTF_BOLD(stderr, "-> If you use this program, please cite the vcfgl paper.\n");
    FPRINTF_BOLD(stderr, "----------------------------------------------------------\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");


    FPRINTF_BOLD_COLOR(stderr, YELLOW, "Option descriptions:\n");
    fprintf(stderr, "     -s, --long-option TYPE [X] _____ Description\n");
    fprintf(stderr, "     -s                               Short option (if any)\n");
    fprintf(stderr, "         --long-option                Long option\n");

    fprintf(stderr, "                       TYPE           Type of the argument value, can be:\n");
    fprintf(stderr, "                                        - INT (integer)\n");
    fprintf(stderr, "                                        - STRING (string)\n");
    fprintf(stderr, "                                        - FILE (filename)\n");
    fprintf(stderr, "                                        - x|y|z (one of the listed values x, y or z)\n");
    fprintf(stderr, "                            [X]       Default argument value (if any)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");

    FPRINTF_BOLD_COLOR(stderr, BLUE, "Usage: setAlleles -i <input> [options]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -h, --help _____________________  Print this help message and exit\n");
    fprintf(stderr, "    -v, --version __________________  Print version and build information and exit\n");
    fprintf(stderr, "\n");
    //fprintf(stderr, "Options:\n");
    FPRINTF_BOLD_COLOR(stderr, BLUE, "Options:\n");
    fprintf(stderr, "    -i, --input FILE ________________ Input BCF file\n");
    fprintf(stderr, "    -a, --alleles FILE ______________ Alleles TSV file\n");
    fprintf(stderr, "    -O, --output-mode [b]|u|z|v _____ Output mode, can be:\n");
    fprintf(stderr, "                                        - b: Compressed BCF (.bcf)\n");
    fprintf(stderr, "                                        - u: uncompressed BCF (.bcf)\n");
    fprintf(stderr, "                                        - z: compressed VCF (.vcf.gz)\n");
    fprintf(stderr, "                                        - v: uncompressed VCF (.vcf)\n");
    fprintf(stderr, "    -o, --output STRING ['output'] __ Output filename prefix\n");
    fprintf(stderr, "\n");
}

args_t* args_get(int argc, char **argv){
    args_t* args = args_init();
    if(argc == 0){
        help_page();
        exit(0);
    }
    while(*argv){
        char* arv = *argv;
        char* val = *(++argv);
        if(strcmp(arv, "-h") == 0 || strcmp(arv, "--help") == 0){
            help_page();
            args_destroy(args);
            exit(0);
        }
        if(strcmp(arv, "-v") == 0 || strcmp(arv, "--version") == 0){
            version_page();
            args_destroy(args);
            exit(0);
        }
        if(strcmp(arv, "-i") == 0 || strcmp(arv, "--input") == 0){
            set_arg_value(&args->input_bcf_file, "--input", val);
        }
        else if(strcmp(arv, "-a") == 0 || strcmp(arv, "--alleles") == 0){
            set_arg_value(&args->alleles_tsv_file, "--alleles", val);
        }
        else if(strcmp(arv, "-O") == 0 || strcmp(arv, "--output-mode") == 0){
            set_arg_value(&args->output_mode, "--output-mode", val);
        }
        else if(strcmp(arv, "-o") == 0 || strcmp(arv, "--output") == 0){
            set_arg_value(&args->output_prefix, "--output", val);
        }  
        ++argv;
    }
    if(args->input_bcf_file == NULL){
        ERROR("Input BCF file is required. Please provide a BCF file using the -i flag (e.g. -i input.bcf)");
    }
    if(args->alleles_tsv_file == NULL){
        ERROR("Alleles TSV file is required. Please provide an alleles TSV file using the -a flag (e.g. -a alleles.tsv)");
    }
    if(args->output_mode == NULL){
        args->output_mode = strdup("b");
    }
    if(args->output_prefix == NULL){
        WARN("Output file prefix is not specified. Setting it to 'output'.\n");
        args->output_prefix = strdup("output");
    }

    switch (*args->output_mode) {
    case 'v':
        args->output_fn = (char*)malloc(strlen(args->output_prefix) + strlen(".vcf") + 1);
        strcpy(args->output_fn, args->output_prefix);
        strcat(args->output_fn, ".vcf");
        free(args->output_mode);
        args->output_mode = strdup("w");
        break;
    case 'b':
        args->output_fn =
            (char*)malloc(strlen(args->output_prefix) + strlen(".bcf") + 1);
        strcpy(args->output_fn, args->output_prefix);
        strcat(args->output_fn, ".bcf");
        free(args->output_mode);
        args->output_mode = strdup("wb");
        break;
    case 'z':
        args->output_fn =
            (char*)malloc(strlen(args->output_prefix) + strlen(".vcf.gz") + 1);
        strcpy(args->output_fn, args->output_prefix);
        strcat(args->output_fn, ".vcf.gz");
        free(args->output_mode);
        args->output_mode = strdup("wz");
        break;
    case 'u':
        args->output_fn =
            (char*)malloc(strlen(args->output_prefix) + strlen(".bcf") + 1);
        strcpy(args->output_fn, args->output_prefix);
        strcat(args->output_fn, ".bcf");
        free(args->output_mode);
        args->output_mode = strdup("wbu");
        break;
    }

    return(args);
}


int main(int argc, char **argv) {
    args_t* args = args_get(--argc, ++argv);
    alleles_t *alleles = load_alleles_tsv(args->alleles_tsv_file);

    // print alleles
    //for (int i = 0; i < alleles->n_sites; ++i) {
    //    printf("%s\t", alleles->ref[i]);
    //    for (int j = 0; j < alleles->n_alts[i]; j++) {
    //        if(j > 0) printf(",");
    //        printf("%s", alleles->alts[i][j]);
    //    }
    //    printf("\n");
    //}

    // input bcf file
    htsFile *in = bcf_open(args->input_bcf_file, "r");
    bcf_hdr_t *in_hdr = bcf_hdr_read(in);
    bcf1_t *in_rec = bcf_init();

    // output bcf file
    htsFile *out = hts_open(args->output_fn, args->output_mode);
    bcf_hdr_t *out_hdr = bcf_hdr_dup(in_hdr);
    bcf1_t *out_rec = bcf_init();
    ASSERT(bcf_hdr_write(out, out_hdr) == 0);

    int n_samples = bcf_hdr_nsamples(in_hdr);

    uint64_t site_idx = 0;
    uint8_t n_new_alleles_at_site;
    uint8_t n_old_alleles_at_site;
    bool ind_is_missing;
    bool found;
    int old_gt_idx;
    int new_gt_idx;
    while (bcf_read(in, in_hdr, in_rec) == 0) {
        ASSERT(site_idx < alleles->n_sites);
        if(site_idx >= alleles->n_sites){
            ERROR("Site index out of bounds. Please make sure the number of lines in the alleles TSV file is the same as the number of sites in the input BCF file.");
        }
        bcf_unpack(in_rec, BCF_UN_ALL);

        // copy input record to output record
        bcf_copy(out_rec, in_rec);
        bcf_unpack(out_rec, BCF_UN_ALL);


        // --------------------------------------------------------------------
        // --> Update alleles <--
        // --------------------------------------------------------------------
        //// print alleles
        //printf("alleles: %s,", alleles->ref[site_idx]);
        //for(uint8_t i = 0; i < alleles->n_alts[site_idx]; ++i){
        //    if(i > 0) printf(",");
        //    printf("%s", alleles->alts[site_idx][i]);
        //}
        //printf("\n");

        // create a comma separated kstring_t of the new alleles
        kstring_t new_alleles = {0, 0, NULL};
        ksprintf(&new_alleles, "%s", alleles->ref[site_idx]);
        kputc(',', &new_alleles);
        for(uint8_t i = 0; i < alleles->n_alts[site_idx]; ++i){
            ksprintf(&new_alleles, "%s", alleles->alts[site_idx][i]);
            if(i < alleles->n_alts[site_idx] - 1){
                kputc(',', &new_alleles);
            }
        }

        bcf_update_alleles_str(out_hdr, out_rec, new_alleles.s);
        ks_free(&new_alleles);
        int32_t ret;

        n_old_alleles_at_site = (uint8_t)in_rec->n_allele;
        n_new_alleles_at_site = alleles->n_alts[site_idx]+1;

        // --------------------------------------------------------------------
        // --> Create old->new and new->old allele idx maps <--
        // --------------------------------------------------------------------
        int olda_to_newa_idx[5] = { -1, -1, -1, -1, -1 };
        //int new_to_old_idx[5] = { -1, -1, -1, -1, -1 };

        //printf("in_rec->n_allele: %d\n", in_rec->n_allele);
        //printf("out_rec->n_allele: %d\n", out_rec->n_allele);

        for(int i = 0; i < in_rec->n_allele; ++i){
            //printf("in_rec->d.allele[%d]: %s\n", i, in_rec->d.allele[i]);
            found=false;
            // match the old allele to the new allele
            for(int j = 0; j < out_rec->n_allele; ++j){
                if(strcmp(in_rec->d.allele[i], alleles->ref[site_idx]) == 0){
                    //printf("old allele %s (idx=%ld) matches new allele %s (idx=0)\n", in_rec->d.allele[i], i, alleles->ref[site_idx]);
                    ASSERT(i < 5);
                    olda_to_newa_idx[i] = 0;
                    //new_to_old_idx[0] = i;
                    found=true;
                }
                else if( j < n_new_alleles_at_site-1 && strcmp(in_rec->d.allele[i], alleles->alts[site_idx][j]) == 0){
                    //printf("old allele %s (idx=%ld) matches new allele %s (idx=%d)\n", in_rec->d.allele[i], i, alleles->alts[site_idx][j], j+1);
                    ASSERT(j+1 < 5);
                    olda_to_newa_idx[i] = j + 1;
                    //new_to_old_idx[j + 1] = i;
                    found=true;
                }
                if(found){
                    break;
                }
            }
        }

        //// print the maps
        //printf("olda_to_newa_idx: ");
        //for(int i = 0; i < n_old_alleles_at_site; ++i){
        //    printf("%d ", olda_to_newa_idx[i]);
        //}
        //printf("\n");

        // --------------------------------------------------------------------
        // --> Create old->new and new->old genotype idx maps <--
        // --------------------------------------------------------------------

        // 0) create oldgt_to_newgt_idx map with max num of gls per sample (15)
        int oldgt_to_newgt_idx[15] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
        int n_old_gts = n_old_alleles_at_site * (n_old_alleles_at_site + 1) / 2;

        //printf("n_old_gts: %d\n", n_old_gts);
        //printf("n_in_gls: %d\n", n_in_gls);

        int n_new_gts = n_new_alleles_at_site * (n_new_alleles_at_site + 1) / 2;
        int n_new_gts_allsamples = ((n_new_alleles_at_site * (n_new_alleles_at_site + 1))/2) * n_samples;

        //printf("n_new_gts: %d\n", n_new_gts);
        //printf("n_new_gts_allsamples: %d\n", n_new_gts_allsamples);

        for(int i=0; i<n_old_gts; ++i){
            int old_a1_idx;
            int old_a2_idx;
            bcf_gt2alleles(i, &old_a1_idx, &old_a2_idx);
            //printf("old genotype idx: %d\n", i);
            //printf("old genotype: %s %s\n", in_rec->d.allele[old_a1_idx], in_rec->d.allele[old_a2_idx]);

            int new_a1_idx = olda_to_newa_idx[old_a1_idx];
            int new_a2_idx = olda_to_newa_idx[old_a2_idx];
            if(new_a1_idx != -1 && new_a2_idx != -1){

                new_gt_idx=bcf_alleles2gt(new_a1_idx, new_a2_idx);
                oldgt_to_newgt_idx[i] = new_gt_idx;

                //printf("new genotype idx: %d\n", new_gt_idx);
                //printf("new genotype: %s %s\n", out_rec->d.allele[new_a1_idx], out_rec->d.allele[new_a2_idx]);
            }
        }

        //printf("\n");
        //printf("\n");
        //// print oldgt_to_newgt_idx
        //printf("oldgt_to_newgt_idx: ");
        //for(int i = 0; i < n_old_gts; ++i){
        //    printf("%d ", oldgt_to_newgt_idx[i]);
        //}
        //printf("\n");
        //printf("\n");

        // --------------------------------------------------------------------
        // --> Update INFO/QS <--
        // --------------------------------------------------------------------

        float *in_info_qs = NULL;
        int n_in_info_qs = 0;
        ret = bcf_get_info_float(in_hdr, in_rec, "QS", &in_info_qs, &n_in_info_qs);
        if(ret > 0){

            //// print old alleles
            //printf("old alleles: ");
            //for(int i = 0; i < in_rec->n_allele; ++i){
            //    printf("%s ", in_rec->d.allele[i]);
            //}
            //printf("\n");

            //// print qs for old alleles
            //printf("qs for old alleles: ");
            //for(int i = 0; i < in_rec->n_allele; ++i){
            //    printf("%f ", in_info_qs[i]);
            //}
            //printf("\n");

            //// print new alleles
            //printf("new alleles: ");
            //for(int i = 0; i < out_rec->n_allele; ++i){
            //    printf("%s ", out_rec->d.allele[i]);
            //}
            //printf("\n");

            // set the new qs using the old qs and the maps
            float *out_info_qs = NULL;
            out_info_qs = (float*)malloc(sizeof(float) * (size_t)n_new_alleles_at_site);
            int32_t n_out_info_qs = n_new_alleles_at_site;
            for(int i = 0; i < n_old_alleles_at_site; ++i){
                if(olda_to_newa_idx[i] != -1){
                    out_info_qs[olda_to_newa_idx[i]] = in_info_qs[i];
                }
            }

            ret = bcf_update_info_float(out_hdr, out_rec, "QS", out_info_qs, n_out_info_qs);

            //// print qs for new alleles
            //printf("qs for new alleles: ");
            //for(int i = 0; i < out_rec->n_allele; ++i){
            //    printf("%f ", out_info_qs[i]);
            //}
            //printf("\n");

            // clean up
            free(out_info_qs);
            out_info_qs = NULL;
        }

        // --------------------------------------------------------------------
        // --> Update FORMAT/GL <--
        // --------------------------------------------------------------------

        // GL (log10-scaled genotype likelihoods)
        float *in_gls = NULL;
        int n_in_gls = 0;

        ret = bcf_get_format_float(in_hdr, in_rec, "GL", &in_gls, &n_in_gls);
        if (ret > 0) {

            float *out_gls = (float*) malloc(sizeof(float) * (size_t)n_new_gts_allsamples);
            ASSERT(out_gls != NULL);
            
            // 1) reconstruct the new gls from the old gls
            for(int s = 0; s < n_samples; ++s){
                for(int i = 0; i < n_old_gts; ++i){
                    new_gt_idx = oldgt_to_newgt_idx[i];
                    if(new_gt_idx != -1){
                        out_gls[s*n_new_gts + new_gt_idx] = in_gls[s*n_old_gts + i];
                        //printf("->\n");
                        //printf("-> s*n_new_gts+i = %d*%d+%d = %d\n", s, n_new_gts, i, s*n_new_gts+i);
                        //printf("-> in_gls[%d*%d+%d] = %f\n", s, n_old_gts, i, in_gls[s*n_old_gts + i]);
                        //printf("-> out_gls[%d*%d+%d] = %f\n", s, n_new_gts, new_gt_idx, out_gls[s*n_new_gts + new_gt_idx]);
                    }
                }
            }

            // 2) renormalize the new gls (N.B. log10-scaled) for each sample
            for(int s = 0; s < n_samples; ++s){
                ind_is_missing = false;
                // first handle missing values (if sample has no genotypes)
                for(int i = s*n_new_gts; i < (s+1)*n_new_gts; ++i){
                    // missing is nan
                    if(isnan(out_gls[i])){
                        bcf_float_set_missing(out_gls[i]);
                        ind_is_missing = true;
                        break;
                    }
                }
                if(ind_is_missing){
                    continue;
                }

                float max_gl = -INFINITY;
                for(int i = s*n_new_gts; i < (s+1)*n_new_gts; ++i){
                    if(out_gls[i] > max_gl){
                        max_gl = out_gls[i];
                    }
                }
                for(int i = s*n_new_gts; i < (s+1)*n_new_gts; ++i){
                    out_gls[i] -= max_gl;
                    //printf("-> out_gls[%d] = %f\n", i, out_gls[i]);
                }
            }

            // print the new gls
            //printf("new gls: ");
            //for(int i = 0; i < n_new_gts_allsamples; ++i){
            //    printf("%f ", out_gls[i]);
            //}
            //printf("\n");

            // print genotype index, genotype, gl
            for(int s=0; s<n_samples; ++s){
                ind_is_missing = false;
                for(int i=0; i<n_old_gts; ++i){
                    // skip if genotype is missing
                    if(bcf_float_is_missing(in_gls[s*n_old_gts + i])){
                        ind_is_missing = true;
                        break;
                    }
                    if(ind_is_missing){
                        continue;
                    }
                    //printf("sample: %d\n", s);
                    //printf("genotype index: %d\n", i);
                    int a1;
                    int a2;
                    bcf_gt2alleles(i, &a1, &a2);
                    //printf("genotype: %s %s\n", in_rec->d.allele[a1], in_rec->d.allele[a2]);
                    //printf("gl: %f\n", out_gls[s*n_new_gts + i]);
                    //printf("%d %d %s%s %f\n", s, i, in_rec->d.allele[a1], in_rec->d.allele[a2], in_gls[s*n_old_gts + i]);
                }
                //printf("\n");
                ind_is_missing = false;
                for(int i=0; i<n_new_gts; ++i){
                    // skip if genotype is missing
                    if(bcf_float_is_missing(out_gls[s*n_new_gts + i])){
                        ind_is_missing = true;
                        break;
                    }
                    if(ind_is_missing){
                        continue;
                    }
                    //printf("sample: %d\n", s);
                    //printf("genotype index: %d\n", i);
                    int a1;
                    int a2;
                    bcf_gt2alleles(i, &a1, &a2);
                    //printf("genotype: %s %s\n", in_rec->d.allele[a1], in_rec->d.allele[a2]);
                    //printf("gl: %f\n", out_gls[s*n_new_gts + i]);
                    //printf("-> s*n_new_gts+i = %d*%d+%d = %d\n", s, n_new_gts, i, s*n_new_gts+i);
                    //printf("%d %d %s%s %f\n", s, i, out_rec->d.allele[a1], out_rec->d.allele[a2], out_gls[s*n_new_gts + i]);
                }
                //printf("\n--------------------------------\n");
            }


            // 3) update the new gls
            ret = bcf_update_format_float(out_hdr, out_rec, "GL", out_gls, n_new_gts_allsamples);

            // clean up
            free(out_gls);
            out_gls = NULL;
        }


        // --------------------------------------------------------------------
        // --> Update FORMAT/PL <--
        // --------------------------------------------------------------------

        // PL (log10-scaled genotype posterior probabilities)
        int32_t *in_pls = NULL;
        int n_in_pls = 0;
        ret = bcf_get_format_int32(in_hdr, in_rec, "PL", &in_pls, &n_in_pls);
        if(ret > 0){
            int32_t *out_pls = (int32_t*) malloc(sizeof(int32_t) * (size_t)n_new_gts_allsamples);
            ASSERT(out_pls != NULL);

            // reconstruct the new pls from the old pls
            for(int s = 0; s < n_samples; ++s){
                for(int i = 0; i < n_old_gts; ++i){
                    new_gt_idx = oldgt_to_newgt_idx[i];
                    if(new_gt_idx != -1){
                        out_pls[s*n_new_gts + new_gt_idx] = in_pls[s*n_old_gts + i];
                    }
                }
            }

            // normalize the new pls
            for(int s = 0; s < n_samples; ++s){
                ind_is_missing = false;
                // first handle missing values (if sample has no genotypes)
                for(int i = s*n_new_gts; i < (s+1)*n_new_gts; ++i){
                    if(out_pls[i] == bcf_int32_missing){
                        out_pls[i] = bcf_int32_missing;
                        ind_is_missing = true;
                        break;
                    }
                }
                if(ind_is_missing){
                    continue;
                }

                float min_pl = INFINITY;
                for(int i = s*n_new_gts; i < (s+1)*n_new_gts; ++i){
                    if(out_pls[i] < min_pl){
                        min_pl = (float)out_pls[i];
                    }
                }
                for(int i = s*n_new_gts; i < (s+1)*n_new_gts; ++i){
                    float temp = (float)out_pls[i] - min_pl;
                    out_pls[i] = (int32_t)temp;
                }
            }

            ret = bcf_update_format_int32(out_hdr, out_rec, "PL", out_pls, n_new_gts_allsamples);

            // clean up
            free(out_pls);
            out_pls = NULL;
        }


        // --------------------------------------------------------------------
        // --> Update FORMAT/GP <--
        // --------------------------------------------------------------------

        // GP (genotype posterior probabilities)
        float *in_gps = NULL;
        int n_in_gps = 0;
        ret = bcf_get_format_float(in_hdr, in_rec, "GP", &in_gps, &n_in_gps);
        if(ret > 0){
            float *out_gps = (float*) malloc(sizeof(float) * (size_t)n_new_gts_allsamples);
            ASSERT(out_gps != NULL);

            // reconstruct the new gps from the old gps
            for(int s = 0; s < n_samples; ++s){
                for(int i = 0; i < n_old_gts; ++i){
                    new_gt_idx = oldgt_to_newgt_idx[i];
                    if(new_gt_idx != -1){
                        out_gps[s*n_new_gts + new_gt_idx] = in_gps[s*n_old_gts + i];
                    }
                }   
            }

            // normalize the new gps
            for(int s = 0; s < n_samples; ++s){
                ind_is_missing = false;
                // first handle missing values (if sample has no genotypes)
                for(int i = s*n_new_gts; i < (s+1)*n_new_gts; ++i){
                    if(isnan(out_gps[i])){
                        bcf_float_set_missing(out_gps[i]);
                        ind_is_missing = true;
                        break;
                    }
                }
                if(ind_is_missing){
                    continue;
                }

                float sum_gps = 0.0;
                for(int i = s*n_new_gts; i < (s+1)*n_new_gts; ++i){
                    sum_gps += out_gps[i];
                }
                for(int i = s*n_new_gts; i < (s+1)*n_new_gts; ++i){
                    out_gps[i] /= sum_gps;
                }
            }

            ret = bcf_update_format_float(out_hdr, out_rec, "GP", out_gps, n_new_gts_allsamples);

            // clean up
            free(out_gps);
            out_gps = NULL;
        }

        // --------------------------------------------------------------------
        // --> Clean up <--
        // --------------------------------------------------------------------

        if (in_gls != NULL) {
            free(in_gls);
            in_gls = NULL;
        }
        if (in_pls != NULL) {
            free(in_pls);
            in_pls = NULL;
        }
        if (in_info_qs != NULL) {
            free(in_info_qs);
            in_info_qs = NULL;
        }
        if (in_gps != NULL) {
            free(in_gps);
            in_gps = NULL;
        }

        ASSERT(bcf_write(out, out_hdr, out_rec) == 0);
        site_idx++;
    }
    ASSERT(site_idx == alleles->n_sites);

    fprintf(stderr, "\n");
    fprintf(stderr, "Program finished successfully!\n");
    fprintf(stderr, "\t-> Processed %ld sites.\n", site_idx);
    fprintf(stderr, "\t-> Output file: %s\n", args->output_fn);


    // clean up
    bcf_destroy(in_rec);
    bcf_hdr_destroy(in_hdr);
    bcf_close(in);

    bcf_destroy(out_rec);
    bcf_hdr_destroy(out_hdr);
    hts_close(out);
    alleles_destroy(alleles);
    args_destroy(args);

    return 0;
}
