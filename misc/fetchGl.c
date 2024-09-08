// v0.1
// isinaltinkaya



#include <stdio.h>
#include <htslib/vcf.h>

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



void usage(void) {
    fprintf(stderr, "Usage: fetchGl -i <input vcf> -gt <genotype to fetch>\n");
    fprintf(stderr, "Example: fetchGl -i input.vcf -gt AA\n");
    return;
}
// program to print gls corresponding to a given gt from a vcf file

int main(int argc, char** argv) {

    char* in_fn = NULL;
    char* gt_to_fetch = NULL;

    --argc;++argv;

    if (0 == argc) {
        usage();
        exit(0);
    }

    while (*argv) {
        char* arv = *argv;
        char* val = *(++argv);
        if ((strcmp("-h", arv) == 0) || (strcmp("--help", arv) == 0)) {
            usage();
            exit(0);
        } else if ((strcmp("-i", arv) == 0) || (strcmp("--input", arv) == 0)) {
            in_fn = strdup(val);
        } else if ((strcmp("-gt", arv) == 0) || (strcmp("--genotype", arv) == 0)) {
            gt_to_fetch = strdup(val);
        } else {
            ERROR("Unknown arg:%s\n", arv);
        }
        ++argv;
    }


    if (in_fn == NULL) {
        ERROR("Input file not specified. Exiting...\n");
    }

    if (gt_to_fetch == NULL) {
        ERROR("Genotype to fetch not specified. Exiting...\n");
    }

    htsFile* vcf = hts_open(in_fn, "r");

    bcf1_t* rec = bcf_init();

    bcf_hdr_t* hdr = bcf_hdr_read(vcf);


    int gtidx = 0;

    const int nSamples = bcf_hdr_nsamples(hdr);

    int n_alleles;

    int gt_to_alleles[2] = { -1, -1 };
    int n_gls = 0;

    while (bcf_read(vcf, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_ALL);

        float* gls = NULL;
        int nvals = 0;
        int size_e = 0;

        n_alleles = rec->n_allele;

        // if (n_alleles < 2) {
        //     fprintf(stderr, "Skipping position %ld since it has less than 2 alleles.\n", rec->pos + 1);
        //     continue;
        // }

        nvals = bcf_get_format_float(hdr, rec, "GL", &gls, &size_e);
        if (nvals <= 0) {
            ERROR("Could not read GL tag. Exiting...\n");
        }
        n_gls = nvals / nSamples;


        char* allele;
        gt_to_alleles[0] = -1;
        gt_to_alleles[1] = -1;

        for (int i = 0;i < n_alleles;++i) {

            allele = rec->d.allele[i];


            if (allele[0] == gt_to_fetch[0]) {
                gt_to_alleles[0] = i;
            }

            if (allele[0] == gt_to_fetch[1]) {
                gt_to_alleles[1] = i;
            }

        }

        if (gt_to_alleles[0] == -1 || gt_to_alleles[1] == -1) {
            fprintf(stderr, "Requested genotype %s not found at position %ld. Skipping...\n", gt_to_fetch, rec->pos + 1);
			fflush(stderr);
            continue;
        }
        fprintf(stderr, "Requested genotype's first allele corresponds to allele %d and second allele corresponds to allele %d at position %ld.\n", gt_to_alleles[0], gt_to_alleles[1], rec->pos + 1);

        gtidx = bcf_alleles2gt(gt_to_alleles[0], gt_to_alleles[1]);



        char glval[50];
        fprintf(stdout, "%ld,", rec->pos + 1);
		fflush(stdout);

        for (int i = 0;i < nSamples;++i) {
            float val = gls[i * n_gls + gtidx];
            if (bcf_float_is_missing(val)) {
                strcpy(glval, "MISSING");
            } else if (bcf_float_is_vector_end(val)) {
                strcpy(glval, "END");
            } else {
                sprintf(glval, "%f", val);
            }

            // check if gl is missing

            // print "Pos,GLs for samples"
            fprintf(stdout, "%s", glval);
            if (i < nSamples - 1) {
                fprintf(stdout, ",");
            } else {
                fprintf(stdout, "\n");
            }
        }



        free(gls);

    }

    free(in_fn);
    free(gt_to_fetch);

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(vcf);

    return(0);

}

