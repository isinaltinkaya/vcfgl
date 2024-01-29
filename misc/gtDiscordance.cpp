// v4
// isinaltinkaya


#include <stdio.h>
#include <htslib/vcf.h>

#include "../shared.h"


// count types
#define COUNT_OVERALL 0
#define HOM_TO_HOM 1
#define HET_TO_HET 2
#define HOM_TO_HET 3
#define HET_TO_HOM 4

#define T_CONCORDANT 0
#define T_DISCORDANT 1

#define MAX_PL 255
#define MIN_PL 0

// #count types
#define N_COUNT_TYPES 5 // 0 1 2 3 4

// #count types that are true
#define N_T_TYPES 3 // 0 1 2 

// (bcftools/mcall.c)
// call->GQs[ismpl] = max <= INT8_MAX ? max : INT8_MAX;
// (stdint.h)
// INT8_MAX 127 
#define MIN_GQ 1
#define MAX_GQ 127
#define gq_arr_input_SIZE 130

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

int alelle_to_int(char* allele) {
	if (allele[0] == 'A' || allele[0] == 'a') {
		return 0;
	} else if (allele[0] == 'C' || allele[0] == 'c') {
		return 1;
	} else if (allele[0] == 'G' || allele[0] == 'g') {
		return 2;
	} else if (allele[0] == 'T' || allele[0] == 't') {
		return 3;
	} else {
		return 4;
	}
}

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
	fprintf(stderr, "  -doGQ 5			per-sample doGQ 3\n");
	fprintf(stderr, "  -doGQ 6			per-sample doGQ 4\n");
	fprintf(stderr, "  -doGQ 7			per-sample doGQ 4, for the sites with missing GQ assume highest GQ value\n");
	fprintf(stderr, "  -doGQ 8			doGQ 7, but for the genotype calling files without GQ tags sets GQ to second lowest PL (assuming lowest PL is 0)\n");
	fprintf(stderr, "  -h, --help       print help\n");
	fprintf(stderr, "\n");

}


int main(int argc, char** argv)
{

	char* in_fn = NULL;
	char* true_fn = NULL;
	char* out_fn = NULL;
	int doGq = 0;
	int perSite = 0;

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
		} else if ((strcasecmp("-perSite", arv) == 0)) {
			perSite = atoi(val);
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

	if (doGq == 0) {
		//
	} else if (doGq == 1) {
		fprintf(stderr, "-doGQ 1. Will print the discordance indicator and GQ scores for each sample at each site.\n");
	} else if (doGq == 2) {
		fprintf(stderr, "-doGQ 2. Will print the discordance indicator, hom/het info and GQ scores for each sample at each site.\n");
	} else if (doGq == 3) {
		fprintf(stderr, "-doGQ 3. Will do gq and attempt to tidy up the output.\n");
	} else if (doGq == 4) {
		fprintf(stderr, "-doGQ 4. Will do gq and attempt to tidy up the output and add true/call hom/het status information.\n");
	} else if (doGq == 5) {
		fprintf(stderr, "-doGQ 5. Will do gq and attempt to tidy up the output for each sample.\n");
	} else if (doGq == 6) {
		fprintf(stderr, "-doGQ 6. Will do gq and attempt to tidy up the output for each sample and add true/call hom/het status information.\n");
	} else if (doGq == 7) {
		fprintf(stderr, "-doGQ 7. Will do gq and attempt to tidy up the output for each sample and add true/call hom/het status information. For the sites with missing GQ assume highest GQ value.\n");
	} else if (doGq == 8) {
		fprintf(stderr, "-doGQ 8. Will do gq and attempt to tidy up the output for each sample and add true/call hom/het status information. For the sites with missing GQ assume highest GQ value. For the genotype calling files without GQ tags sets GQ to second lowest PL (assuming lowest PL is 0).\n");
	} else {
		ERROR("Unknown doGq:%d\n", doGq);
	}

	fprintf(stderr, "N.B. Samples are assumed to appear in the same order and quantity in both files.\n");

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

	int32_t ngt1 = 0, * gt_arr_true = NULL, ngt_arr_true = 0;
	int32_t ngt2 = 0, * gt_arr_input = NULL, ngt_arr_input = 0;

	int32_t ngq = 0, * gq_arr_input = NULL, ngq_arr_input = 0;

	int32_t npl = 0, * pl_arr_input = NULL, npl_arr_input = 0;

	int gqScore = -1;


	int* nSitesDiscordant = (int*)calloc(nSamples, sizeof(int));

	// number of sites (per individual) where the call file has the site but ind i is missing
	int* nSites_callmis = (int*)calloc(nSamples, sizeof(int));

	int* nSites_compared_forSample = (int*)calloc(nSamples, sizeof(int));
	int trueb1 = 0, trueb2 = 0, inputb1 = 0, inputb2 = 0;

	double missingness_rate = 0.0;
	double discordance_rate = 0.0;
	double concordance_rate = 0.0;
	int nSitesRetained = 0;

	// assume both files have one same contig
	ASSERT(0 == strcmp(bcf_hdr_id2name(hdrTrue, recTrue->rid), bcf_hdr_id2name(hdrInput, recInput->rid)));

	int nSites1 = 0;
	int nSites2 = 0;
	int nSitesfTrueNotInfInput = 0;

	int dcType = 0;
	int refToCallType = 0;

	const int nSamplesAndTotal = nSamples + 1;
	const int total_across_samples = nSamples; // index to get the total across all samples

	// \def gqCounts [i][j][k][m]  = number of sites with REF_TO_CALL type j with gq score k for sample i with discordant/concordant info m
	// e.g. gqCounts[0][HOM_TO_HET][30][T_DISCORDANT] = number of sites with gq score 30 for sample 0(first sample) that is truth:hom call:het and type discordant
	// e.g. gqCounts[total_across_samples][HOM_TO_HET][30][T_DISCORDANT] = number of sites with gq score 30 for all samples in total that is truth:hom call:het and type discordant
	int**** gqCounts = NULL;
	gqCounts = (int****)malloc(nSamplesAndTotal * sizeof(int***)); // last element is for sample_i=TOTAL FOR ALL SAMPLES
	ASSERT(gqCounts != NULL);

	for (int i = 0;i < nSamplesAndTotal;++i) {
		gqCounts[i] = (int***)malloc(N_COUNT_TYPES * sizeof(int**));
		ASSERT(gqCounts[i] != NULL);
		for (int j = 0;j < N_COUNT_TYPES;++j) {
			gqCounts[i][j] = (int**)malloc(gq_arr_input_SIZE * sizeof(int*));
			ASSERT(gqCounts[i][j] != NULL);
			for (int k = 0;k < gq_arr_input_SIZE;++k) {
				gqCounts[i][j][k] = (int*)malloc(2 * sizeof(int));
				ASSERT(gqCounts[i][j][k] != NULL);
			}
		}
	}
	for (int i = 0;i < nSamplesAndTotal;++i) {
		for (int j = 0;j < N_COUNT_TYPES;++j) {
			for (int k = 0;k < gq_arr_input_SIZE;++k) {
				gqCounts[i][j][k][T_DISCORDANT] = 0;
				gqCounts[i][j][k][T_CONCORDANT] = 0;
			}
		}
	}

	int ret1 = -1;
	int ret2 = -1;


	ret1 = bcf_read(fTrue, hdrTrue, recTrue);
	nSites1++;

	bool hasGq = (bcf_hdr_id2int(hdrInput, BCF_DT_ID, "GQ") >= 0);
	if (doGq == 7) {

		ASSERT(hasGq == true);
	}

	bool hasPl = (bcf_hdr_id2int(hdrInput, BCF_DT_ID, "PL") >= 0);
	if ((doGq == 8) && (!hasGq)) {
		if (!hasPl) {
			ERROR("-doGQ 8 is only supported for genotype calling files with GQ or PL tags");
		}
	}

	bool hasSymbolicAlt = false;



	while (0 == (ret2 = bcf_read(fInput, hdrInput, recInput))) {

		nSites2++;

		while (recInput->pos > recTrue->pos) {


			ret1 = bcf_read(fTrue, hdrTrue, recTrue);
			nSites1++;

			ASSERT(ret1 == 0); // expected since true file should always contain the call file sites

		}


		ASSERT(recInput->pos == recTrue->pos);

		// prep for new rec
		ngt1 = 0;
		ngt_arr_true = 0;
		free(gt_arr_true);
		gt_arr_true = NULL;
		ngt2 = 0;
		ngt_arr_input = 0;
		free(gt_arr_input);
		gt_arr_input = NULL;
		ngq = 0;
		ngq_arr_input = 0;
		if (gq_arr_input != NULL) {
			free(gq_arr_input);
			gq_arr_input = NULL;
		}


		// begin new rec



		ASSERT(0 == bcf_unpack(recTrue, BCF_UN_STR));
		ASSERT(0 == bcf_unpack(recInput, BCF_UN_STR));

		ngt1 = bcf_get_genotypes(hdrTrue, recTrue, &gt_arr_true, &ngt_arr_true);
		ASSERT(ngt1 > 0);
		ngt2 = bcf_get_genotypes(hdrInput, recInput, &gt_arr_input, &ngt_arr_input);
		ASSERT(ngt2 > 0);

		// get genotype quality scores from the input/call file
		if (doGq > 0) {
			ngq = bcf_get_format_int32(hdrInput, recInput, "GQ", &gq_arr_input, &ngq_arr_input);

			if (doGq == 7) {

				if (gq_arr_input == NULL) {
					// if GQ is missing for site, then it is a site identified as invar, assume highest GQ value
					gq_arr_input = (int32_t*)malloc(nSamples * sizeof(int32_t));
					for (int i = 0; i < nSamples; i++)
					{
						gq_arr_input[i] = MAX_GQ;
					}
				}

			} else if (doGq == 8) {


				if (hasGq) {

					if (gq_arr_input == NULL) {
						// if GQ is missing for site, then it is a site identified as invar, assume highest GQ value
						gq_arr_input = (int32_t*)malloc(nSamples * sizeof(int32_t));
						for (int i = 0; i < nSamples; i++)
						{
							// missing case is already handled since we check for missing gt later in program before using the gq_arr_input
							// so no need to check and handle it here
							gq_arr_input[i] = MAX_GQ;
						}
					}

				} else {


					ASSERT(gq_arr_input == NULL);


					// if ind is nonmissing and npl = 1 and/or pl is missing, assuma highest GQ value

					// if GQ is missing entirely, set GQ to second lowest PL (assuming lowest PL is 0)
					// if second lowest PL > max GQ, set to max GQ

					gq_arr_input = (int32_t*)malloc(nSamples * sizeof(int32_t));
					for (int i = 0; i < nSamples; i++)
					{
						gq_arr_input[i] = 0;
					}



					int n_alleles = recInput->n_allele;
					hasSymbolicAlt = false;
					for (int i = 1; i < n_alleles; ++i)
					{
						if (recInput->d.allele[i][0] == '<' && recInput->d.allele[i][1] == '*' && recInput->d.allele[i][2] == '>') {
							hasSymbolicAlt = true;
							break;
						} else if (recInput->d.allele[i][0] == '<' && recInput->d.allele[i][1] == 'N' && recInput->d.allele[i][2] == 'O' && recInput->d.allele[i][3] == 'N' && recInput->d.allele[i][4] == '_' && recInput->d.allele[i][5] == 'R' && recInput->d.allele[i][6] == 'E' && recInput->d.allele[i][7] == 'F' && recInput->d.allele[i][8] == '>') {
							hasSymbolicAlt = true;
							break;
						}
					}


					// find second lowest PL for each sample
					npl = bcf_get_format_int32(hdrInput, recInput, "PL", &pl_arr_input, &npl_arr_input);
					ASSERT(npl > 0);

					int nGenotypes = 0;

					do {

						if (n_alleles == 0) {
							NEVER;
						} else if (n_alleles == 1) {
							ASSERT(hasSymbolicAlt == false);
							nGenotypes = 1;
							for (int i = 0; i < nSamples; ++i)
							{
								gq_arr_input[i] = MAX_GQ;
							}
							break;
						} else if (n_alleles == 2) {
							nGenotypes = 3;
						} else if (n_alleles == 3) {
							nGenotypes = 6;
						} else if (n_alleles == 4) {
							nGenotypes = 10;
						} else if (n_alleles == 5) {
							nGenotypes = 15;
						} else {
							NEVER;
						}

						int secondbestpl = MAX_PL;
						bool foundZero = false;

						for (int i = 0; i < nSamples; ++i)
						{
							foundZero = false;
							secondbestpl = MAX_PL;
							// fprintf(stderr, "\n-> nGenotypes=%d", nGenotypes);


							for (int g = 0;g < nGenotypes;++g) {
								int32_t plval = pl_arr_input[i * nGenotypes + g];

								if (bcf_int32_missing == plval) {
									// missing case is already handled since we check for missing gt later in program before using the gq_arr_input
									// so no need to check and handle it here
									foundZero = true;
									secondbestpl = MAX_PL;
									// continue;




								} else if (0 == plval) {
									// if (foundZero) {
										// ERROR("Found more than one zero PL for sample %d at site %d", i, nSites2);
									// }
									foundZero = true;
								} else if (bcf_int32_vector_end == plval) {
									NEVER;
								} else if (pl_arr_input[i * nGenotypes + g] < secondbestpl) {
									secondbestpl = pl_arr_input[i * nGenotypes + g];
								}


							}

							gq_arr_input[i] = secondbestpl > MAX_GQ ? MAX_GQ : secondbestpl;
							ASSERT(foundZero == true);

						}


					} while (0);





				}


			} else {
				ASSERT(ngq > 0);
				ASSERT(ngq == nSamples);

			}
		}


		for (int i = 0; i < nSamples; i++)
		{



			int32_t* ptr_true = gt_arr_true + i * 2;
			int32_t* ptr_input = gt_arr_input + i * 2;

			if (ptr_true[0] == bcf_int32_vector_end || ptr_true[1] == bcf_int32_vector_end) {
				NEVER;
			}
			if (ptr_input[0] == bcf_int32_vector_end || ptr_input[1] == bcf_int32_vector_end) {
				NEVER;
			}


			if (bcf_gt_is_missing(ptr_true[0]) || bcf_gt_is_missing(ptr_true[1])) {
				NEVER;
			}


			if (bcf_gt_is_missing(ptr_input[0]) || bcf_gt_is_missing(ptr_input[1])) {
				nSites_callmis[i]++;
				continue;
			}






			dcType = -1;
			nSites_compared_forSample[i]++;

			trueb1 = refToInt[(int)*recTrue->d.allele[bcf_gt_allele(ptr_true[0])]];
			trueb2 = refToInt[(int)*recTrue->d.allele[bcf_gt_allele(ptr_true[1])]];

			inputb1 = refToInt[(int)*recInput->d.allele[bcf_gt_allele(ptr_input[0])]];
			inputb2 = refToInt[(int)*recInput->d.allele[bcf_gt_allele(ptr_input[1])]];


			ASSERT(trueb1 >= 0);
			ASSERT(trueb1 <= 3);
			ASSERT(trueb2 >= 0);
			ASSERT(trueb2 <= 3);

			ASSERT(inputb1 >= 0);
			ASSERT(inputb1 <= 3);
			ASSERT(inputb2 >= 0);
			ASSERT(inputb2 <= 3);


			refToCallType = -1;

			if (trueb1 == trueb2) {
				// HOM_TO_
				if (inputb1 == inputb2) {
					// HOM_TO_HOM
					refToCallType = HOM_TO_HOM;
					if (trueb1 == inputb1) {
						dcType = T_CONCORDANT;
					} else {
						dcType = T_DISCORDANT;
					}
				} else {
					// HOM_TO_HET
					refToCallType = HOM_TO_HET;
					dcType = T_DISCORDANT;
				}
			} else {
				// HET_TO_
				if (inputb1 == inputb2) {
					// HET_TO_HOM
					refToCallType = HET_TO_HOM;
					dcType = T_DISCORDANT;
				} else {
					// HET_TO_HET
					refToCallType = HET_TO_HET;
					if (trueb1 == inputb1) {
						if (trueb2 == inputb2) {
							dcType = T_CONCORDANT;
						} else {
							dcType = T_DISCORDANT;
						}
					} else {
						if (trueb2 == inputb2) {
							dcType = T_DISCORDANT;
						} else {
							if (trueb1 == inputb2) {
								if (trueb2 == inputb1) {
									dcType = T_CONCORDANT;
								} else {
									dcType = T_DISCORDANT;
								}
							} else {
								dcType = T_DISCORDANT;
							}
						}
					}
				}
			}

			if (NULL != gq_arr_input) {
				gqScore = gq_arr_input[i];
				ASSERT(gqScore <= MAX_GQ);
				ASSERT(gqScore > 0); // can be 0 iff no call, and we already skip nocalls before this
			}

			if (T_DISCORDANT == dcType) {
				nSitesDiscordant[i]++;
			}

			if (doGq > 0) {
				gqCounts[i][COUNT_OVERALL][gqScore][dcType]++;
				gqCounts[i][refToCallType][gqScore][dcType]++;
				gqCounts[total_across_samples][COUNT_OVERALL][gqScore][dcType]++;
				gqCounts[total_across_samples][refToCallType][gqScore][dcType]++;
			} else {
				gqCounts[i][COUNT_OVERALL][0][dcType]++;
				gqCounts[i][refToCallType][0][dcType]++;
				gqCounts[total_across_samples][COUNT_OVERALL][0][dcType]++;
				gqCounts[total_across_samples][refToCallType][0][dcType]++;
			}

			if (2 == doGq) {
				// colnames(d)<-c("Site","Ind","isDiscordant","refToCallType","GQ")
				fprintf(out_fp, "%d\t%d\t%d\t%d\t%d\n", nSites1, i, dcType, refToCallType, gqScore);
			} else if (1 == doGq) {
				// colnames(d)<-c("Site","Ind","isDiscordant","GQ")
				fprintf(out_fp, "%d\t%d\t%d\t%d\n", nSites1, i, dcType, gqScore);
			}

			if (1 == perSite) {
				// colnames(d)<-c("Pos","Ind","isDiscordant","refToCallType")
				fprintf(stdout, "%ld\t%d\t%d\t%d\n", recTrue->pos + 1, i, dcType, refToCallType);
			}

		} // samples loop

	}
	while (0 == (ret1 = bcf_read(fTrue, hdrTrue, recTrue))) {
		nSites1++;
	}

	nSitesfTrueNotInfInput = nSites1 - nSites2;

	fprintf(stderr, "Comparing the genotypes in the truth file %s and the call file %s\n", true_fn, in_fn);
	fprintf(stderr, "Number of sites in truth file: %d\n", nSites1);
	fprintf(stderr, "Number of sites in the truth file not in the call file: %d\n", nSitesfTrueNotInfInput);
	fprintf(stderr, "Number of sites used in comparison: %d\n", nSites2);

	nSitesRetained = nSites1 - nSitesfTrueNotInfInput;
	ASSERT(nSitesRetained == nSites2);

	if (3 == doGq) {
		// k = 1 since GQ==0 is empty; 0 GQ is already skipped due to missing GT
		for (int k = 1;k < gq_arr_input_SIZE;++k) {
			// colnames(d)<-c("GQ","nDiscordant","nConcordant")
			fprintf(out_fp, "%d\t%d\t%d\n", k,
				gqCounts[total_across_samples][COUNT_OVERALL][k][T_DISCORDANT],
				gqCounts[total_across_samples][COUNT_OVERALL][k][T_CONCORDANT]);
		}

	} else if (5 == doGq) {
		// per-sample doGq 3
		for (int i = 0; i < nSamples; i++)
		{
			for (int k = 1;k < gq_arr_input_SIZE;++k) {
				// N.B. nSitesComparedForSample == nSitesInTotal (in truth file including sites with no non0 dp) - nSitesDP0forSample
				// colnames(d)<-c("Sample","GQ","nDiscordant","nConcordant","nSitesComparedForSample")
				// N.B. nSitesComparedForSample is the overall nSitesComparedForSample. So same value is repeated.
				// perGQ_nSitesComparedForSample can be obtained using nDiscordant + nConcordant 
				fprintf(out_fp, "%d\t%d\t%d\t%d\t%d\n",
					i,
					k,
					gqCounts[i][COUNT_OVERALL][k][T_DISCORDANT],
					gqCounts[i][COUNT_OVERALL][k][T_CONCORDANT],
					nSites_compared_forSample[i]
				);
			}
		}

	} else if (4 == doGq) {

		for (int k = 1;k < gq_arr_input_SIZE;++k) {
			fprintf(out_fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
				k,
				gqCounts[total_across_samples][COUNT_OVERALL][k][T_DISCORDANT],
				gqCounts[total_across_samples][HOM_TO_HOM][k][T_DISCORDANT],
				gqCounts[total_across_samples][HOM_TO_HET][k][T_DISCORDANT],
				gqCounts[total_across_samples][HET_TO_HOM][k][T_DISCORDANT],
				gqCounts[total_across_samples][HET_TO_HET][k][T_DISCORDANT],
				gqCounts[total_across_samples][COUNT_OVERALL][k][T_CONCORDANT],
				gqCounts[total_across_samples][HOM_TO_HOM][k][T_CONCORDANT],
				gqCounts[total_across_samples][HET_TO_HET][k][T_CONCORDANT]);
		}

	} else if (6 == doGq) {

		// per-sample doGq 4
		for (int i = 0; i < nSamples; i++)
		{
			// colnames(d)<-c("Sample","GQ","nDiscordant","nHomToHomDiscordant","nHomToHetDiscordant","nHetToHomDiscordant","nHetToHetDiscordant","nConcordant","nHomToHomConcordant","nHetToHetConcordant","nSitesComparedForSample")
			// N.B. nOverall == nSitesComparedForSample == nSitesInTotal - nSitesDP0forSample
			for (int k = 1;k < gq_arr_input_SIZE;++k) {
				fprintf(out_fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
					i,
					k,
					gqCounts[i][COUNT_OVERALL][k][T_DISCORDANT],
					gqCounts[i][HOM_TO_HOM][k][T_DISCORDANT],
					gqCounts[i][HOM_TO_HET][k][T_DISCORDANT],
					gqCounts[i][HET_TO_HOM][k][T_DISCORDANT],
					gqCounts[i][HET_TO_HET][k][T_DISCORDANT],
					gqCounts[i][COUNT_OVERALL][k][T_CONCORDANT],
					gqCounts[i][HOM_TO_HOM][k][T_CONCORDANT],
					gqCounts[i][HET_TO_HET][k][T_CONCORDANT],
					nSites_compared_forSample[i]
				);
			}
		}

	} else if (7 == doGq) {

		// per-sample doGq 4 with missing GQ assumed to be highest GQ
		for (int i = 0; i < nSamples; i++)
		{
			// colnames(d)<-c("Sample","GQ","nDiscordant","nHomToHomDiscordant","nHomToHetDiscordant","nHetToHomDiscordant","nHetToHetDiscordant","nConcordant","nHomToHomConcordant","nHetToHetConcordant","nSitesComparedForSample")
			// N.B. nOverall == nSitesComparedForSample == nSitesInTotal - nSitesDP0forSample
			for (int k = 1;k < gq_arr_input_SIZE;++k) {
				fprintf(out_fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
					i,
					k,
					gqCounts[i][COUNT_OVERALL][k][T_DISCORDANT],
					gqCounts[i][HOM_TO_HOM][k][T_DISCORDANT],
					gqCounts[i][HOM_TO_HET][k][T_DISCORDANT],
					gqCounts[i][HET_TO_HOM][k][T_DISCORDANT],
					gqCounts[i][HET_TO_HET][k][T_DISCORDANT],
					gqCounts[i][COUNT_OVERALL][k][T_CONCORDANT],
					gqCounts[i][HOM_TO_HOM][k][T_CONCORDANT],
					gqCounts[i][HET_TO_HET][k][T_CONCORDANT],
					nSites_compared_forSample[i]
				);
			}
		}
	} else if (8 == doGq) {


		for (int i = 0; i < nSamples; i++)
		{
			// colnames(d)<-c("Sample","GQ","nDiscordant","nHomToHomDiscordant","nHomToHetDiscordant","nHetToHomDiscordant","nHetToHetDiscordant","nConcordant","nHomToHomConcordant","nHetToHetConcordant","nSitesComparedForSample")
			// N.B. nOverall == nSitesComparedForSample == nSitesInTotal - nSitesDP0forSample
			for (int k = 1;k < gq_arr_input_SIZE;++k) {
				fprintf(out_fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
					i,
					k,
					gqCounts[i][COUNT_OVERALL][k][T_DISCORDANT],
					gqCounts[i][HOM_TO_HOM][k][T_DISCORDANT],
					gqCounts[i][HOM_TO_HET][k][T_DISCORDANT],
					gqCounts[i][HET_TO_HOM][k][T_DISCORDANT],
					gqCounts[i][HET_TO_HET][k][T_DISCORDANT],
					gqCounts[i][COUNT_OVERALL][k][T_CONCORDANT],
					gqCounts[i][HOM_TO_HOM][k][T_CONCORDANT],
					gqCounts[i][HET_TO_HET][k][T_CONCORDANT],
					nSites_compared_forSample[i]
				);
			}
		}

	} else if (0 == doGq) {
		for (int i = 0; i < nSamples; i++)
		{

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

			missingness_rate = (double)(1.0 - ((double)nSites_compared_forSample[i] / (double)nSites1));
			discordance_rate = (double)nSitesDiscordant[i] / (double)nSites_compared_forSample[i];
			concordance_rate = (double)(nSites_compared_forSample[i] - nSitesDiscordant[i]) / (double)nSites_compared_forSample[i];

			T_Hom = gqCounts[i][HOM_TO_HOM][0][T_CONCORDANT];
			T_Het = gqCounts[i][HET_TO_HET][0][T_CONCORDANT];

			F_HomToHom = gqCounts[i][HOM_TO_HOM][0][T_DISCORDANT];
			F_HomToHet = gqCounts[i][HOM_TO_HET][0][T_DISCORDANT];
			F_HetToHom = gqCounts[i][HET_TO_HOM][0][T_DISCORDANT];
			F_HetToHet = gqCounts[i][HET_TO_HET][0][T_DISCORDANT];

			T_Hom_rate = (double)T_Hom / (double)nSites_compared_forSample[i];
			T_Het_rate = (double)T_Het / (double)nSites_compared_forSample[i];
			F_HomToHom_rate = (double)F_HomToHom / (double)nSites_compared_forSample[i];
			F_HomToHet_rate = (double)F_HomToHet / (double)nSites_compared_forSample[i];
			F_HetToHom_rate = (double)F_HetToHom / (double)nSites_compared_forSample[i];
			F_HetToHet_rate = (double)F_HetToHet / (double)nSites_compared_forSample[i];

			fprintf(out_fp, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n",
				hdrTrue->samples[i],
				nSites1,
				nSitesRetained,
				nSites_compared_forSample[i],
				nSites_callmis[i],
				nSitesDiscordant[i],
				nSitesfTrueNotInfInput,
				nSites_compared_forSample[i] - nSitesDiscordant[i],
				missingness_rate,
				discordance_rate,
				concordance_rate,
				T_Hom,
				T_Het,
				F_HomToHom,
				F_HomToHet,
				F_HetToHom,
				F_HetToHet,
				T_Hom_rate,
				T_Het_rate,
				F_HomToHom_rate,
				F_HomToHet_rate,
				F_HetToHom_rate,
				F_HetToHet_rate);

			// header:
			// "Sample","nSitesTotal","nSitesRetained","nSitesCompared","nSitesCallMis","nDiscordantSites","nConcordantSites","nSitesinTrueNotCall","MissingnessRate","DiscordanceRate","ConcordanceRate","T_Hom","T_Het","F_HomToHom","F_HomToHet","F_HetToHom","F_HetToHet","T_Hom_rate","T_Het_rate","F_HomToHom_rate","F_HomToHet_rate","F_HetToHom_rate","F_HetToHet_rate"

			// sanity check
			// nSites_compared_forSample is the number of sites with nonmissing data for sample i
			ASSERT(nSites_compared_forSample[i] == nSitesRetained - nSites_callmis[i]);


		}

	}


	for (int i = 0;i < nSamplesAndTotal;++i) {
		for (int j = 0;j < N_COUNT_TYPES;++j) {
			for (int k = 0;k < gq_arr_input_SIZE;++k) {
				free(gqCounts[i][j][k]);
				gqCounts[i][j][k] = NULL;
			}
			free(gqCounts[i][j]);
			gqCounts[i][j] = NULL;
		}
		free(gqCounts[i]);
		gqCounts[i] = NULL;
	}
	free(gqCounts);
	gqCounts = NULL;

	free(nSitesDiscordant);
	nSitesDiscordant = NULL;
	free(nSites_callmis);
	nSites_callmis = NULL;
	free(nSites_compared_forSample);
	nSites_compared_forSample = NULL;
	free(gt_arr_true);
	gt_arr_true = NULL;
	free(gt_arr_input);
	gt_arr_input = NULL;
	if (gq_arr_input != NULL) {
		free(gq_arr_input);
		gq_arr_input = NULL;
	}
	if (pl_arr_input != NULL) {
		free(pl_arr_input);
		pl_arr_input = NULL;
	}
	bcf_destroy(recTrue);
	bcf_destroy(recInput);
	bcf_hdr_destroy(hdrTrue);
	bcf_hdr_destroy(hdrInput);
	bcf_close(fTrue);
	bcf_close(fInput);

	fprintf(stderr, "\n\nFinished running gtDiscordance\n");
	fprintf(stderr, "\t-> Output file: %s\n", out_fn);
	fclose(out_fp);
	free(in_fn);
	in_fn = NULL;
	free(true_fn);
	true_fn = NULL;
	free(out_fn);
	out_fn = NULL;


	return 0;
}
