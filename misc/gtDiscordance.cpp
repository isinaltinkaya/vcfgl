// v2
// isinaltinkaya


#include <stdio.h>
#include <htslib/vcf.h>

#include "../shared.h"

#define COUNT_OVERALL 0
#define HOM_TO_HOM 1
#define HET_TO_HET 2
#define HOM_TO_HET 3
#define HET_TO_HOM 4
#define N_COUNT_TYPES 5 // 0 1 2 3 4
#define N_T_TYPES 3 // 0 1 2

// (bcftools/mcall.c)
// call->GQs[ismpl] = max <= INT8_MAX ? max : INT8_MAX;
// (stdint.h)
// INT8_MAX 127 
#define MAX_GQ 127
#define GQ_ARR_SIZE 130

// maps a->0,A->0,c->1,C->1,g->2,G->2,t->3,T=>3,n->4,N->5
int refToInt[256] = {
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 15
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 31
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 47
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 63
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 79
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 95
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 111
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 127
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 143
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 159
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 175
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 191
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 207
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 223
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 239
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4  // 255
};

void usage(void) {

    fprintf(stderr, "Usage: ./gtDiscordance -t <truth.bcf> -i <call.bcf> -o <output.tsv>\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -t, --truth <file>		truth file\n");
    fprintf(stderr, "  -i, --input <file>		input file, to be compared with truth file, typically a genotype call file\n");
    fprintf(stderr, "  -o, --output <file>		output\n");
    fprintf(stderr, "  -doGQ 1			print the discordance indicator and GQ scores for each sample at each site\n");
    fprintf(stderr, "  -doGQ 2			print the discordance indicator and GQ scores for each sample at each site, with extra info\n");
    fprintf(stderr, "  -doGQ 3			print tidy gq summary\n");
    fprintf(stderr, "  -doGQ 4			print tidy gq summary with hom/het details\n");
    fprintf(stderr, "  -h, --help			print help\n");
    fprintf(stderr, "\n");


}

int main(int argc, char** argv)
{

    char* in_fn = NULL;
    char* true_fn = NULL;
    char* out_fn = NULL;
    int doGq = 0;

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
        } else if ((strcmp("-t", arv) == 0) || (strcmp("--truth", arv) == 0)) {
            true_fn = strdup(val);
        } else if ((strcmp("-o", arv) == 0) || (strcmp("--output", arv) == 0)) {
            out_fn = strdup(val);
        } else if ((strcasecmp("-doGQ", arv) == 0)) {
            doGq = atoi(val);
        } else {
            ERROR("Unknown arg:%s\n", arv);
        }
        ++argv;
    }

    //check number of args

    if (NULL == in_fn) {
        ERROR("Input file not specified. Please use -i/--input <file>.\n");
    }
    if (NULL == true_fn) {
        ERROR("Truth file not specified. Please use -t/--truth <file>.\n");
    }
    if (NULL == out_fn) {
        ERROR("Output file not specified. Please use -o/--output <file>.\n");
    }
    fprintf(stderr, "%s", out_fn);

    if (doGq == 1) {
        fprintf(stderr, "-doGQ 1. Will print the discordance indicator and GQ scores for each sample at each site.\n");
    } else if (doGq == 2) {
        fprintf(stderr, "-doGQ 2. Will print the discordance indicator, hom/het info and GQ scores for each sample at each site.\n");
    } else if (doGq == 3) {
        fprintf(stderr, "-doGQ 3. Will do gq and attempt to tidy up the output.\n");
    } else if (doGq == 4) {
        fprintf(stderr, "-doGQ 4. Will do gq and attempt to tidy up the output and add true/call hom/het status information.\n");
    }

    htsFile* fTrue = bcf_open(true_fn, "r");
    ASSERT(fTrue != NULL);

    htsFile* fInput = bcf_open(in_fn, "r");
    ASSERT(fInput != NULL);

    // output file
    FILE* out_fp = NULL;
    out_fp = fopen(out_fn, "w");
    if (out_fp == NULL) {
        ERROR("Could not open output file %s for writing", out_fn);
    }

    bcf_hdr_t* hdrTrue = bcf_hdr_read(fTrue);
    bcf_hdr_t* hdrInput = bcf_hdr_read(fInput);
    bcf1_t* recTrue = bcf_init();
    bcf1_t* recInput = bcf_init();

    int nSamples = bcf_hdr_nsamples(hdrTrue);

    int32_t ngt1 = 0, * gt_arr1 = NULL, ngt_arr1 = 0;
    int32_t ngt2 = 0, * gt_arr2 = NULL, ngt_arr2 = 0;

    int32_t ngq = 0, * gq_arr = NULL, ngq_arr = 0;


    int* nSitesDiscordant = (int*)calloc(nSamples, sizeof(int));

    // number of sites (per individual) where the call file has the site but ind i is missing
    int* nSites_callmis = (int*)calloc(nSamples, sizeof(int));

    int* nSites_compared_forSample = (int*)calloc(nSamples, sizeof(int));
    int f1b1 = 0, f1b2 = 0, f2b1 = 0, f2b2 = 0;
    int f1a1 = 0, f1a2 = 0, f2a1 = 0, f2a2 = 0;
    int sidx1 = 0, sidx2 = 0;

    double missingness_rate = 0.0;
    double discordance_rate = 0.0;
    double concordance_rate = 0.0;
    int nSitesRetained = 0;

    // assume both files have one same contig
    ASSERT(0 == strcmp(bcf_hdr_id2name(hdrTrue, recTrue->rid), bcf_hdr_id2name(hdrInput, recInput->rid)));

    int nSites1 = 0;
    int nSites2 = 0;
    int nSitesfTrueNotInfInput = 0;

    // call: correct, truegt=callgt: hom
    // output id=0
    int T_Hom = 0;

    // call: correct, truegt=callgt: het
    // output id=1
    int T_Het = 0;

    // call: incorrect, truegt: hom, callgt: hom
    // e.g. truegt: A/A, callgt: C/C
    // output id=2
    int F_HomToHom = 0;

    // call: incorrect, truegt: hom, callgt: het
    // output id=3
    int F_HomToHet = 0;

    // call: incorrect, truegt: het, callgt: hom
    // output id=4
    int F_HetToHom = 0;

    // call: incorrect, truegt: het, callgt: het
    // e.g. truegt: A/C, callgt: A/T
    // N.B. unphased therefore A/C == C/A
    // output id=5
    int F_HetToHet = 0;

    double T_Hom_rate = 0.0;
    double T_Het_rate = 0.0;
    double F_HomToHom_rate = 0.0;
    double F_HomToHet_rate = 0.0;
    double F_HetToHom_rate = 0.0;
    double F_HetToHet_rate = 0.0;

    int f1_gt_binary_1 = 0;
    int f1_gt_binary_2 = 0;
    int f2_gt_binary_1 = 0;
    int f2_gt_binary_2 = 0;

    // 0: hom, 1: het
    int f1_homhet_state = 0;
    int f2_homhet_state = 0;

    int is_discordant = 0;


    // \def gqDiscordantCounts[i][j] = number of sites type of i with gq score j
    // e.g. gqDiscordantCounts[HOM_TO_HET][10] = number of discordant sites with truth:hom call:het and gq score 10
    // e.g. gqConcordant[0][30] = number of concordant sites with gq score 30
    int gqDiscordantCounts[N_COUNT_TYPES][GQ_ARR_SIZE] = { { 0 } };
    int gqConcordantCounts[N_T_TYPES][GQ_ARR_SIZE] = { { 0 } };


    while (bcf_read(fTrue, hdrTrue, recTrue) == 0)
    {
        nSites1++;

        if (bcf_read(fInput, hdrInput, recInput) != 0)
        {

            // fInput finished
            if (recInput == NULL)
            {
                // continue with next rec in fTrue
                continue;
            } else
            {
                nSitesfTrueNotInfInput++;
                // continue with next rec in fTrue
                continue;
            }
            NEVER;
        }
        nSites2++;

        ASSERT(0 == bcf_unpack(recTrue, BCF_UN_STR));
        ASSERT(0 == bcf_unpack(recInput, BCF_UN_STR));


        ngt1 = bcf_get_genotypes(hdrTrue, recTrue, &gt_arr1, &ngt_arr1);
        ngt2 = bcf_get_genotypes(hdrInput, recInput, &gt_arr2, &ngt_arr2);

        for (int i = 0; i < nSamples; i++)
        {
            sidx1 = 2 * i;
            sidx2 = sidx1 + 1;

            f1_gt_binary_1 = gt_arr1[sidx1];
            f1_gt_binary_2 = gt_arr1[sidx2];
            f2_gt_binary_1 = gt_arr2[sidx1];
            f2_gt_binary_2 = gt_arr2[sidx2];

            (f1_gt_binary_1 == f1_gt_binary_2) ? (f1_homhet_state = 0) : (f1_homhet_state = 1);
            (f2_gt_binary_1 == f2_gt_binary_2) ? (f2_homhet_state = 0) : (f2_homhet_state = 1);


            f1a1 = bcf_gt_allele(f1_gt_binary_1);
            f1a2 = bcf_gt_allele(f1_gt_binary_2);

            f2a1 = bcf_gt_allele(f2_gt_binary_1);
            f2a2 = bcf_gt_allele(f2_gt_binary_2);

            // true gt file (fTrue) can never have missing
            ASSERT(f1a1 != -1);
            ASSERT(f1a2 != -1);

            if (-1 == f2a1)
            {
                if (-1 == f2a2)
                {
                    nSites_callmis[i]++;
                    continue;
                } else
                {
                    NEVER;
                }
            }

            f1b1 = refToInt[*recTrue->d.allele[bcf_gt_allele(gt_arr1[sidx1])]];
            f1b2 = refToInt[*recTrue->d.allele[bcf_gt_allele(gt_arr1[sidx2])]];

            f2b1 = refToInt[*recInput->d.allele[bcf_gt_allele(gt_arr2[sidx1])]];
            f2b2 = refToInt[*recInput->d.allele[bcf_gt_allele(gt_arr2[sidx2])]];

            ASSERT(f1b1 >= 0);
            ASSERT(f1b1 <= 3);
            ASSERT(f1b2 >= 0);
            ASSERT(f1b2 <= 3);

            ASSERT(f2b1 >= 0);
            ASSERT(f2b1 <= 3);
            ASSERT(f2b2 >= 0);
            ASSERT(f2b2 <= 3);

            // get genotype quality scores from the input/call file
            if (doGq > 0) {
                ngq = bcf_get_format_int32(hdrInput, recInput, "GQ", &gq_arr, &ngq_arr);
                if (ngq <= 0) {
                    ERROR("No GQ scores for sample %d\n", i);
                }
            }

            nSites_compared_forSample[i]++;

            is_discordant = 0;

            if (f1b1 != f2b1)
            {
                if (f1b1 != f2b2)
                {
                    is_discordant = 1;
                } else
                {
                    // f2b2 used;exclude from check
                    // now check for f1b2
                    if (f1b2 != f2b1)
                    {
                        is_discordant = 1;
                    }
                }
            } else
            {
                // f2b1 used;exclude from check
                // now check for f1b2
                if (f1b2 != f2b2)
                {
                    is_discordant = 1;
                }
            }



            ASSERT(gq_arr[i] < 128);
            if (1 == is_discordant) {
                if (3 == doGq || 4 == doGq) {
                    gqDiscordantCounts[COUNT_OVERALL][gq_arr[i]]++;
                }
                if (0 == f1_homhet_state) {
                    if (0 == f2_homhet_state) {
                        F_HomToHom++;
                        if (2 == doGq) {
                            fprintf(out_fp, "%d\t%d\t%d\t2\t%d\n", nSites1, i, is_discordant, gq_arr[i]);
                        }
                        if (3 == doGq || 4 == doGq) {
                            gqDiscordantCounts[HOM_TO_HOM][gq_arr[i]]++;
                        }
                    } else {
                        F_HomToHet++;
                        if (2 == doGq) {
                            fprintf(out_fp, "%d\t%d\t%d\t3\t%d\n", nSites1, i, is_discordant, gq_arr[i]);
                        }
                        if (3 == doGq || 4 == doGq) {
                            gqDiscordantCounts[HOM_TO_HET][gq_arr[i]]++;
                        }
                    }
                } else {
                    if (0 == f2_homhet_state) {
                        F_HetToHom++;
                        if (2 == doGq) {
                            fprintf(out_fp, "%d\t%d\t%d\t4\t%d\n", nSites1, i, is_discordant, gq_arr[i]);
                        }
                        if (3 == doGq || 4 == doGq) {
                            gqDiscordantCounts[HET_TO_HOM][gq_arr[i]]++;
                        }
                    } else {
                        F_HetToHet++;
                        if (2 == doGq) {
                            fprintf(out_fp, "%d\t%d\t%d\t5\t%d\n", nSites1, i, is_discordant, gq_arr[i]);
                        }
                        if (3 == doGq || 4 == doGq) {
                            gqDiscordantCounts[HET_TO_HET][gq_arr[i]]++;
                        }
                    }
                }
                nSitesDiscordant[i]++;
            } else {
                if (3 == doGq) {
                    gqConcordantCounts[COUNT_OVERALL][gq_arr[i]]++;
                }
                // concordant
                if (0 == f1_homhet_state) {
                    T_Hom++;
                    if (2 == doGq) {
                        fprintf(out_fp, "%d\t%d\t%d\t0\t%d\n", nSites1, i, is_discordant, gq_arr[i]);
                    }
                    if (3 == doGq || 4 == doGq) {
                        gqConcordantCounts[HOM_TO_HOM][gq_arr[i]]++;
                    }
                } else {
                    if (2 == doGq) {
                        fprintf(out_fp, "%d\t%d\t%d\t1\t%d\n", nSites1, i, is_discordant, gq_arr[i]);
                    }
                    T_Het++;
                    if (3 == doGq || 4 == doGq) {
                        gqConcordantCounts[HET_TO_HET][gq_arr[i]]++;
                    }
                }
            }

            // colnames(d)<-c("Site","Ind","isDiscordant","GQ")
            if (1 == doGq) {
                fprintf(out_fp, "%d\t%d\t%d\t%d\n", nSites1, i, is_discordant, gq_arr[i]);
                // }else if (2==doGq){
                    // colnames(d)<-c("Site","Ind","isDiscordant","fromTo","GQ")
                    // fprintf(out_fp, "%d\n", gq_arr[i]);
            }

        } // samples loop
    }

    fprintf(stderr, "Comparing the genotypes in the truth file %s and the call file %s\n", true_fn, in_fn);
    fprintf(stderr, "Number of sites in truth file: %d\n", nSites1);
    fprintf(stderr, "Number of sites in the truth file not in the call file: %d\n", nSitesfTrueNotInfInput);
    fprintf(stderr, "Number of sites used in comparison: %d\n", nSites2);

    nSitesRetained = nSites1 - nSitesfTrueNotInfInput;
    ASSERT(nSitesRetained == nSites2);

    if (3 == doGq) {
        for (int i = 0;i < GQ_ARR_SIZE;++i) {
            // colnames(d)<-c("GQ","nDiscordant","nConcordant")
            fprintf(out_fp, "%d\t%d\t%d\n", i, gqDiscordantCounts[COUNT_OVERALL][i], gqConcordantCounts[COUNT_OVERALL][i]);
        }
    } else if (4 == doGq) {
        for (int i = 0;i < GQ_ARR_SIZE;++i) {
            fprintf(out_fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", i, gqDiscordantCounts[COUNT_OVERALL][i], gqDiscordantCounts[HOM_TO_HOM][i], gqDiscordantCounts[HOM_TO_HET][i], gqDiscordantCounts[HET_TO_HOM][i], gqDiscordantCounts[HET_TO_HET][i], gqConcordantCounts[COUNT_OVERALL][i], gqConcordantCounts[HOM_TO_HOM][i], gqConcordantCounts[HET_TO_HET][i]);
        }

    } else {
        for (int i = 0; i < nSamples; i++)
        {

            missingness_rate = (double)(1.0 - ((double)nSites_compared_forSample[i] / (double)nSites1));

            discordance_rate = (double)nSitesDiscordant[i] / (double)nSites_compared_forSample[i];

            concordance_rate = (double)(nSites_compared_forSample[i] - nSitesDiscordant[i]) / (double)nSites_compared_forSample[i];

            T_Hom_rate = (double)T_Hom / (double)nSites_compared_forSample[i];
            T_Het_rate = (double)T_Het / (double)nSites_compared_forSample[i];
            F_HomToHom_rate = (double)F_HomToHom / (double)nSites_compared_forSample[i];
            F_HomToHet_rate = (double)F_HomToHet / (double)nSites_compared_forSample[i];
            F_HetToHom_rate = (double)F_HetToHom / (double)nSites_compared_forSample[i];
            F_HetToHet_rate = (double)F_HetToHet / (double)nSites_compared_forSample[i];

            fprintf(out_fp, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", hdrTrue->samples[i], nSites1, nSitesRetained, nSites_compared_forSample[i], nSites_callmis[i], nSitesDiscordant[i], nSites_compared_forSample[i] - nSitesDiscordant[i], missingness_rate, discordance_rate, concordance_rate, T_Hom, T_Het, F_HomToHom, F_HomToHet, F_HetToHom, F_HetToHet, T_Hom_rate, T_Het_rate, F_HomToHom_rate, F_HomToHet_rate, F_HetToHom_rate, F_HetToHet_rate);

            // header:
            // "Sample","nSitesTotal","nSitesRetained","nSitesCompared","nSitesCallMis","nDiscordantSites","nNonDiscordantSites","MissingnessRate","DiscordanceRate","ConcordanceRate","T_Hom","T_Het","F_HomToHom","F_HomToHet","F_HetToHom","F_HetToHet","T_Hom_rate","T_Het_rate","F_HomToHom_rate","F_HomToHet_rate","F_HetToHom_rate","F_HetToHet_rate"

            // sanity check
            // nSites_compared_forSample is the number of sites with nonmissing data for sample i
            ASSERT(nSites_compared_forSample[i] == nSitesRetained - nSites_callmis[i]);
        }

    }


    free(nSitesDiscordant);
    nSitesDiscordant = NULL;
    free(nSites_callmis);
    nSites_callmis = NULL;
    free(nSites_compared_forSample);
    nSites_compared_forSample = NULL;
    free(gt_arr1);
    gt_arr1 = NULL;
    free(gt_arr2);
    gt_arr2 = NULL;
    if (gq_arr != NULL) {
        free(gq_arr);
        gq_arr = NULL;
    }
    bcf_destroy(recTrue);
    bcf_destroy(recInput);
    bcf_hdr_destroy(hdrTrue);
    bcf_hdr_destroy(hdrInput);
    bcf_close(fTrue);
    bcf_close(fInput);

    fclose(out_fp);
    free(in_fn);
    in_fn = NULL;
    free(true_fn);
    true_fn = NULL;
    free(out_fn);
    out_fn = NULL;


    return 0;
}
