#include "bcf_utils.h"
#include "io.h"



inline float* bcf_tag_alloc_max_size(enum bcf_tag t, const double init_val, simRecord* sim) {
    const int size = sim->max_size_bcf_tag_number[bcf_tags[t].n];
    float* arr = (float*)malloc(size * sizeof(float));
    for (int i = 0; i < size; ++i) {
        arr[i] = init_val;
    }
    ASSERT(NULL != arr);
    return(arr);
}


inline float* bcf_tag_alloc_max_size(enum bcf_tag t, const float init_val, simRecord* sim) {
    const int size = sim->max_size_bcf_tag_number[bcf_tags[t].n];
    float* arr = (float*)malloc(size * sizeof(float));
    for (int i = 0; i < size; ++i) {
        arr[i] = init_val;
    }
    ASSERT(NULL != arr);
    return (arr);
}

inline int32_t* bcf_tag_alloc_max_size(enum bcf_tag t, const int32_t init_val, simRecord* sim) {
    const int size = sim->max_size_bcf_tag_number[bcf_tags[t].n];
    int32_t* arr = (int32_t*)malloc(size * sizeof(int32_t));
    for (int i = 0; i < size; ++i) {
        arr[i] = init_val;
    }
    ASSERT(NULL != arr);
    return (arr);
}


union {
    uint32_t i = bcf_float_missing;
    float f;
} bcf_float_missing_union = { .i = bcf_float_missing };
float bcf_float_missing_union_f = bcf_float_missing_union.f;


// sets simRecord->hdr and simRecord->truth_hdr
simRecord::simRecord(bcf_hdr_t* in_hdr) {

    this->_nBasesPerSample = BUFSIZE_NBASES;

    this->nSamples = bcf_hdr_nsamples(in_hdr);
    this->nContigs = in_hdr->n[BCF_DT_CTG];

    if (this->nContigs != 1) {
        ERROR("Only one contig simulation is supported at the moment. If this is a feature you need please open an issue on GitHub.");
    }

    this->nHaplotypes = this->nSamples * SIM_PLOIDY;

    if (args->printPileup) {
        this->pileup = (kstring_t*)malloc(sizeof(kstring_t));
        this->pileup->s = NULL;
        this->pileup->l = 0;
        this->pileup->m = 0;
    }


    if (0 == args->preCalc) {
        this->base_qScores = (int**)malloc(this->nSamples * sizeof(double*));
        for (int s = 0;s < this->nSamples;++s) {
            this->base_qScores[s] = (int*)malloc(this->_nBasesPerSample * sizeof(int));
            for (int i = 0;i < this->_nBasesPerSample;++i) {
                this->base_qScores[s][i] = -1;
            }
        }

        if (PROGRAM_WILL_ADJUST_QS) {
            this->adj_base_qScores = (int**)malloc(this->nSamples * sizeof(double*));
            for (int s = 0;s < this->nSamples;++s) {
                this->adj_base_qScores[s] = (int*)malloc(this->_nBasesPerSample * sizeof(int));
                for (int i = 0;i < this->_nBasesPerSample;++i) {
                    this->adj_base_qScores[s][i] = -1;
                }
            }
        }
    }

    if ((0 == args->preCalc) && (1 == args->usePreciseGlError)) {
        this->base_error_probs = (double**)malloc(this->nSamples * sizeof(double*));
        for (int s = 0;s < this->nSamples;++s) {
            this->base_error_probs[s] = (double*)malloc(this->_nBasesPerSample * sizeof(double));
            for (int i = 0;i < this->_nBasesPerSample;++i) {
                this->base_error_probs[s][i] = -1.0;
            }
        }
    }

    this->bases = (int**)malloc(this->nSamples * sizeof(int*));
    for (int s = 0;s < this->nSamples;++s) {
        this->bases[s] = (int*)malloc(this->_nBasesPerSample * sizeof(int));
        for (int i = 0;i < this->_nBasesPerSample;++i) {
            this->bases[s][i] = -1;
        }
    }


    this->rec = bcf_init();
    if (args->doGVCF) {
        ASSERT(NULL != (this->gvcfd = gvcfData_init()));
    }

    // ACGT: A, C, G, T, <*> (non-ref)
    // assumption: MAX_NALLELES == 5
    this->acgt2alleles = (int*)malloc(MAX_NALLELES * sizeof(int));
    this->acgt2alleles[0] = -1;
    this->acgt2alleles[1] = -1;
    this->acgt2alleles[2] = -1;
    this->acgt2alleles[3] = -1;
    this->acgt2alleles[4] = -1;


    this->alleles2acgt = (int*)malloc(MAX_NALLELES * sizeof(int));
    this->alleles2acgt[0] = -1;
    this->alleles2acgt[1] = -1;
    this->alleles2acgt[2] = -1;
    this->alleles2acgt[3] = -1;
    this->alleles2acgt[4] = -1;


    // acgt arrays always only contain A,C,G,T (and no non-ref i.e. <*> or <NON_REF>)
    // therefore always 4
    this->acgt_info_ad_arr = (int32_t*)malloc(4 * sizeof(int32_t));
    this->acgt_info_ad_arr[0] = 0;
    this->acgt_info_ad_arr[1] = 0;
    this->acgt_info_ad_arr[2] = 0;
    this->acgt_info_ad_arr[3] = 0;

    this->acgt_fmt_ad_arr = (int32_t*)malloc(4 * this->nSamples * sizeof(int32_t));
    this->acgt_fmt_adf_arr = (int32_t*)malloc(4 * this->nSamples * sizeof(int32_t));
    this->acgt_fmt_adr_arr = (int32_t*)malloc(4 * this->nSamples * sizeof(int32_t));
    this->acgt_fmt_qsum_arr = (int32_t*)malloc(4 * this->nSamples * sizeof(int32_t));
    this->acgt_fmt_qsum_sq_arr = (int32_t*)malloc(4 * this->nSamples * sizeof(int32_t));
    for (int i = 0; i < nSamples * 4; ++i) {
        this->acgt_fmt_ad_arr[i] = 0;
        this->acgt_fmt_adf_arr[i] = 0;
        this->acgt_fmt_adr_arr[i] = 0;
        this->acgt_fmt_qsum_arr[i] = 0;
        this->acgt_fmt_qsum_sq_arr[i] = 0;
    }

    // A_SIM_FORWARD_STRAND,A_SIM_REVERSE_STRAND, C_SIM_FORWARD_STRAND, C_SIM_REVERSE_STRAND, G_SIM_FORWARD_STRAND, G_SIM_REVERSE_STRAND, T_SIM_FORWARD_STRAND, T_SIM_REVERSE_STRAND
    this->acgt_n_bases_forI16 = (int*)malloc(8 * sizeof(int));
    this->acgt_n_bases_forI16[0] = 0; // A_SIM_FORWARD_STRAND
    this->acgt_n_bases_forI16[1] = 0; // A_SIM_REVERSE_STRAND
    this->acgt_n_bases_forI16[2] = 0; // C_SIM_FORWARD_STRAND
    this->acgt_n_bases_forI16[3] = 0; // C_SIM_REVERSE_STRAND
    this->acgt_n_bases_forI16[4] = 0; // G_SIM_FORWARD_STRAND
    this->acgt_n_bases_forI16[5] = 0; // G_SIM_REVERSE_STRAND
    this->acgt_n_bases_forI16[6] = 0; // T_SIM_FORWARD_STRAND
    this->acgt_n_bases_forI16[7] = 0; // T_SIM_REVERSE_STRAND


    this->acgt_sum_taildist = (float*)malloc(4 * sizeof(float));
    this->acgt_sum_taildist[0] = 0.0;
    this->acgt_sum_taildist[1] = 0.0;
    this->acgt_sum_taildist[2] = 0.0;
    this->acgt_sum_taildist[3] = 0.0;

    this->acgt_sum_taildist_sq = (float*)malloc(4 * sizeof(float));
    this->acgt_sum_taildist_sq[0] = 0.0;
    this->acgt_sum_taildist_sq[1] = 0.0;
    this->acgt_sum_taildist_sq[2] = 0.0;
    this->acgt_sum_taildist_sq[3] = 0.0;

    // set the maximum size of each bcf tag array
    // then init all current tag array sizes to max
    // if trimming not used then these are not modified so they remain equal to
    // max
    this->max_size_bcf_tag_number = (int*)malloc(N_ENUM_BCF_TAG_NUMBER * sizeof(int));
    this->current_size_bcf_tag_number = (int*)malloc(N_ENUM_BCF_TAG_NUMBER * sizeof(int));

    this->max_size_bcf_tag_number[FMT_NUMBER_1] = this->nSamples;
    this->current_size_bcf_tag_number[FMT_NUMBER_1] = this->max_size_bcf_tag_number[FMT_NUMBER_1];

    this->max_size_bcf_tag_number[FMT_NUMBER_GT] = this->nHaplotypes;
    this->current_size_bcf_tag_number[FMT_NUMBER_GT] = this->max_size_bcf_tag_number[FMT_NUMBER_GT];

    this->max_size_bcf_tag_number[FMT_NUMBER_G] = this->nSamples * MAX_NGTS;
    // variable size depending on the number of possible genotypes
    this->current_size_bcf_tag_number[FMT_NUMBER_G] = -1;

    this->max_size_bcf_tag_number[FMT_NUMBER_R] = this->nSamples * (MAX_NALLELES - 1);
    // variable size depending on the number of possible alleles at each site
    this->current_size_bcf_tag_number[FMT_NUMBER_R] = -1;

    this->max_size_bcf_tag_number[FMT_NUMBER_R_WITH_NONREF] = this->nSamples * MAX_NALLELES;
    // variable size depending on the number of possible alleles at each site
    this->current_size_bcf_tag_number[FMT_NUMBER_R_WITH_NONREF] = -1;

    this->max_size_bcf_tag_number[INFO_NUMBER_1] = 1;
    this->current_size_bcf_tag_number[INFO_NUMBER_1] = 1;

    this->max_size_bcf_tag_number[INFO_NUMBER_G] = MAX_NGTS;
    this->current_size_bcf_tag_number[INFO_NUMBER_G] = -1;

    this->max_size_bcf_tag_number[INFO_NUMBER_R] = (MAX_NALLELES - 1);
    // variable size depending on the number of possible alleles at each site
    this->current_size_bcf_tag_number[INFO_NUMBER_R] = -1;

    this->max_size_bcf_tag_number[INFO_NUMBER_R_WITH_NONREF] = MAX_NALLELES;
    // variable size depending on the number of possible alleles at each site
    this->current_size_bcf_tag_number[INFO_NUMBER_R_WITH_NONREF] = -1;

    this->max_size_bcf_tag_number[INFO_NUMBER_16] = 16;
    this->current_size_bcf_tag_number[INFO_NUMBER_16] = 16;

    this->set_hdr(in_hdr);

    if (args->mps_depths != NULL) {
        if (args->n_mps_depths != this->nSamples) {
            ERROR("Number of depths provided in the depths file (%d) does not match the number of samples (%d) in the input VCF file.", args->n_mps_depths, this->nSamples);
        }
    }

    this->gt_arr = NULL;  // allocated in-place


    // arr array: used for writing to tag
    // always created, only added if addTYPE==1

    this->gl_arr = bcf_tag_alloc_max_size(GL, -0.0, this);
    this->fmt_dp_arr = bcf_tag_alloc_max_size(FMT_DP, 0, this);
    this->info_dp_arr = bcf_tag_alloc_max_size(INFO_DP, 0, this);

    if (args->addGP) {
        this->gp_arr = bcf_tag_alloc_max_size(GP, MINGP, this);
    }

    if (args->addPL) {
        this->pl_arr = bcf_tag_alloc_max_size(PL, MAXPL, this);
    }

    if (args->addQS) {
        this->qs_arr = bcf_tag_alloc_max_size(QS, 0.0, this);
    }

    if (args->addI16) {
        this->i16_arr = bcf_tag_alloc_max_size(I16, 0.0, this);
    }

    if (args->addFormatAD || args->addI16) {
        this->fmt_ad_arr = bcf_tag_alloc_max_size(FMT_AD, 0, this);
    }

    if (args->addFormatADF) {
        this->fmt_adf_arr = bcf_tag_alloc_max_size(FMT_ADF, 0, this);
    }

    if (args->addFormatADR) {
        this->fmt_adr_arr = bcf_tag_alloc_max_size(FMT_ADR, 0, this);
    }

    if (args->addInfoAD) {
        this->info_ad_arr = bcf_tag_alloc_max_size(INFO_AD, 0, this);
    }

    if (args->addInfoADF) {
        this->info_adf_arr = bcf_tag_alloc_max_size(INFO_ADF, 0, this);
    }

    if (args->addInfoADR) {
        this->info_adr_arr = bcf_tag_alloc_max_size(INFO_ADR, 0, this);
    }

}

simRecord::~simRecord() {


    ASSERT(args != NULL);

    for (int s = 0;s < this->nSamples;++s) {
        if (NULL != this->base_error_probs) {
            free(this->base_error_probs[s]);
            this->base_error_probs[s] = NULL;
        }
        if (NULL != this->base_qScores) {
            free(this->base_qScores[s]);
            this->base_qScores[s] = NULL;
        }
        if (NULL != this->adj_base_qScores) {
            free(this->adj_base_qScores[s]);
            this->adj_base_qScores[s] = NULL;
        }
        free(this->bases[s]);
        this->bases[s] = NULL;
    }
    if (NULL != this->base_error_probs) {
        free(this->base_error_probs);
        this->base_error_probs = NULL;
    }
    if (NULL != this->base_qScores) {
        free(this->base_qScores);
        this->base_qScores = NULL;
    }
    if (NULL != this->adj_base_qScores) {
        free(this->adj_base_qScores);
        this->adj_base_qScores = NULL;
    }
    free(this->bases);
    this->bases = NULL;


    if (args->printPileup) {
        free(this->pileup->s);
        this->pileup->s = NULL;
        free(this->pileup);
        this->pileup = NULL;
    }

    free(this->acgt_fmt_ad_arr);
    this->acgt_fmt_ad_arr = NULL;
    free(this->acgt_fmt_adf_arr);
    this->acgt_fmt_adf_arr = NULL;
    free(this->acgt_fmt_adr_arr);
    this->acgt_fmt_adr_arr = NULL;
    free(this->acgt_fmt_qsum_arr);
    this->acgt_fmt_qsum_arr = NULL;
    free(this->acgt_fmt_qsum_sq_arr);
    this->acgt_fmt_qsum_sq_arr = NULL;

    this->alleles.m = 0;
    this->alleles.l = 0;
    free(this->alleles.s);
    this->alleles.s = NULL;

    free(this->acgt_info_ad_arr);
    this->acgt_info_ad_arr = NULL;

    free(this->acgt_n_bases_forI16);
    this->acgt_n_bases_forI16 = NULL;

    free(this->acgt_sum_taildist);
    this->acgt_sum_taildist = NULL;

    free(this->acgt_sum_taildist_sq);
    this->acgt_sum_taildist_sq = NULL;

    bcf_hdr_destroy(this->hdr);

    if (NULL != this->truth_hdr) {
        bcf_hdr_destroy(this->truth_hdr);
    }

    free(this->gt_arr);
    this->gt_arr = NULL;

    free(this->gl_arr);
    this->gl_arr = NULL;

    free(this->fmt_dp_arr);
    this->fmt_dp_arr = NULL;

    free(this->info_dp_arr);
    this->info_dp_arr = NULL;

    if (NULL != this->info_ad_arr) {
        free(this->info_ad_arr);
        this->info_ad_arr = NULL;
    }

    if (NULL != this->info_adf_arr) {
        free(this->info_adf_arr);
        this->info_adf_arr = NULL;
    }

    if (NULL != this->info_adr_arr) {
        free(this->info_adr_arr);
        this->info_adr_arr = NULL;
    }

    if (NULL != this->gp_arr) {
        free(this->gp_arr);
        this->gp_arr = NULL;
    }

    if (NULL != this->pl_arr) {
        free(this->pl_arr);
        this->pl_arr = NULL;
    }

    if (NULL != this->qs_arr) {
        free(this->qs_arr);
        this->qs_arr = NULL;
    }

    if (NULL != this->i16_arr) {
        free(this->i16_arr);
        this->i16_arr = NULL;
    }

    if (NULL != this->fmt_ad_arr) {
        free(this->fmt_ad_arr);
        this->fmt_ad_arr = NULL;
    }

    if (NULL != this->fmt_adf_arr) {
        free(this->fmt_adf_arr);
        this->fmt_adf_arr = NULL;
    }

    if (NULL != this->fmt_adr_arr) {
        free(this->fmt_adr_arr);
        this->fmt_adr_arr = NULL;
    }

    free(this->max_size_bcf_tag_number);
    this->max_size_bcf_tag_number = NULL;

    free(this->current_size_bcf_tag_number);
    this->current_size_bcf_tag_number = NULL;

    free(this->alleles2acgt);
    this->alleles2acgt = NULL;
    free(this->acgt2alleles);
    this->acgt2alleles = NULL;

    if (this->gvcfd != NULL) {
        gvcfData_destroy(this->gvcfd);
    }

}




void simRecord::add_tags() {

    if (args->addFormatDP) {
        ASSERT(0 == (bcf_update_format_int32(this->hdr, this->rec, "DP", this->fmt_dp_arr, this->current_size_bcf_tag_number[bcf_tags[FMT_DP].n])));
    }

    if (args->addInfoDP) {
        ASSERT(0 ==
            (bcf_update_info_int32(
                this->hdr, this->rec, "DP", this->info_dp_arr,
                this->current_size_bcf_tag_number[bcf_tags[INFO_DP].n])));
    }


    // FMT_NUMBER_G with NON_REF (if set)
    if (args->addGL) {
        ASSERT(0 == (bcf_update_format_float(this->hdr, this->rec, "GL", this->gl_arr, this->current_size_bcf_tag_number[bcf_tags[GL].n])));
    }

    if (args->addPL) {
        ASSERT(0 == (bcf_update_format_int32(
            this->hdr, this->rec, "PL", this->pl_arr,
            this->current_size_bcf_tag_number[bcf_tags[PL].n])));
    }

    if (args->addGP) {
        ASSERT(0 == (bcf_update_format_float(
            this->hdr, this->rec, "GP", this->gp_arr,
            this->current_size_bcf_tag_number[bcf_tags[GP].n])));
    }

    // FMT_NUMBER_G_NO_NONREF
    if (args->addQS) {
        ASSERT(0 == (bcf_update_info_float(
            this->hdr, this->rec, "QS", this->qs_arr,
            this->current_size_bcf_tag_number[bcf_tags[QS].n])));
    }

    if (args->addI16) {
        ASSERT(0 == (bcf_update_info_float(this->hdr, this->rec, "I16",
            this->i16_arr, 16)));
    }

    if (args->addFormatAD) {
        ASSERT(0 == (bcf_update_format_int32(this->hdr, this->rec, "AD", this->fmt_ad_arr, this->current_size_bcf_tag_number[bcf_tags[FMT_AD].n])));
    }

    if (args->addFormatADF) {
        ASSERT(0 ==
            (bcf_update_format_int32(
                this->hdr, this->rec, "ADF", this->fmt_adf_arr,
                this->current_size_bcf_tag_number[bcf_tags[FMT_ADF].n])));
    }

    if (args->addFormatADR) {
        ASSERT(0 ==
            (bcf_update_format_int32(
                this->hdr, this->rec, "ADR", this->fmt_adr_arr,
                this->current_size_bcf_tag_number[bcf_tags[FMT_ADR].n])));
    }

    if (args->addInfoAD) {
        ASSERT(0 ==
            (bcf_update_info_int32(
                this->hdr, this->rec, "AD", this->info_ad_arr,
                this->current_size_bcf_tag_number[bcf_tags[INFO_AD].n])));
    }

    if (args->addInfoADF) {
        ASSERT(0 ==
            (bcf_update_info_int32(
                this->hdr, this->rec, "ADF", this->info_adf_arr,
                this->current_size_bcf_tag_number[bcf_tags[INFO_ADF].n])));
    }

    if (args->addInfoADR) {
        (bcf_update_info_int32(
            this->hdr, this->rec, "ADR", this->info_adr_arr,
            this->current_size_bcf_tag_number[bcf_tags[INFO_ADR].n]));
    }

}



void simRecord::set_hdr(bcf_hdr_t* in_hdr) {

    this->hdr = bcf_hdr_dup(in_hdr);

    char* DATE_TAG = NULL;
    ASSERT(asprintf(&DATE_TAG, "##fileDate=%s", args->datetime) > 0);
    ASSERT(0 == bcf_hdr_append(this->hdr, DATE_TAG));
    free(DATE_TAG);
    DATE_TAG = NULL;

    char* SOURCE_VERSION_TAG;
    ASSERT(asprintf(&SOURCE_VERSION_TAG, "##source=%s", args->versionInfo) > 0);
    ASSERT(0 == bcf_hdr_append(this->hdr, SOURCE_VERSION_TAG));
    free(SOURCE_VERSION_TAG);
    SOURCE_VERSION_TAG = NULL;

    if (args->addI16) {
        char* SOURCE_I16_WARNING;
        ASSERT(asprintf(&SOURCE_I16_WARNING, "##source=[vcfgl] I16 tag simulation is still under development. Please use with caution and report any issues") > 0);
        ASSERT(0 == bcf_hdr_append(this->hdr, SOURCE_I16_WARNING));
        free(SOURCE_I16_WARNING);
        SOURCE_I16_WARNING = NULL;
    }

    char* SOURCE_TAG = NULL;
    ASSERT(asprintf(&SOURCE_TAG, "##source=%s", args->command) > 0);
    ASSERT(0 == bcf_hdr_append(this->hdr, SOURCE_TAG));
    free(SOURCE_TAG);
    SOURCE_TAG = NULL;

    if (args->printTruth) {
        this->truth_hdr = bcf_hdr_dup(this->hdr);
    }
    bcf_hdr_remove(this->hdr, BCF_HL_FMT, "GT");

    if (PROGRAM_WILL_ADD_STAR) {
        ASSERT(0 == bcf_hdr_append(this->hdr, "##ALT=<ID=*,Description=\"Symbolic alternate allele representing any possible alternative allele at this location\">"));

    } else if (PROGRAM_WILL_ADD_NONREF) {
        ASSERT(0 == bcf_hdr_append(this->hdr, "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">"));
    }

    if (1 == args->doGVCF) {


        bcf_hdr_append(this->hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of this non-variant block\">\n");

        // - bcftools 1.18:
        // bcf_hdr_append(this->hdr, "##INFO=<ID=MinDP,Number=1,Type=Integer,Description=\"Minimum per-sample depth in this gVCF block\">");
        // - bcftools 1.18-31-g4014f7e2-dirty: 
        bcf_hdr_append(this->hdr, "##INFO=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP observed within the GVCF block\">");

    }




    if (args->addFormatDP) {
        ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[FMT_DP].hdr));
    }

    if (args->addInfoDP) {
        ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[INFO_DP].hdr));
    }

    if (args->addGL) {
        ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[GL].hdr));
    }

    if (args->addPL) {
        ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[PL].hdr));
    }

    if (args->addGP) {
        ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[GP].hdr));
    }

    if (args->addQS) {
        ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[QS].hdr));
    }

    if (args->addI16) {
        ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[I16].hdr));
    }

    if (args->addFormatAD) {
        ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[FMT_AD].hdr));
    }
    if (args->addFormatADF) {
        ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[FMT_ADF].hdr));
    }
    if (args->addFormatADR) {
        ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[FMT_ADR].hdr));
    }

    if (args->addInfoAD) {
        ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[INFO_AD].hdr));
    }
    if (args->addInfoADF) {
        ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[INFO_ADF].hdr));
    }
    if (args->addInfoADR) {
        ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[INFO_ADR].hdr));
    }
}


void simRecord::expand_arrays(const int new_size) {

    int old_nBasesPerSample = this->_nBasesPerSample;
    this->_nBasesPerSample = new_size;
    if (NULL != this->base_error_probs) {
        for (int s = 0; s < nSamples; s++) {
            this->base_error_probs[s] = (double*)realloc(this->base_error_probs[s], this->_nBasesPerSample * sizeof(double));
            for (int i = old_nBasesPerSample;i < this->_nBasesPerSample;++i) {
                this->base_error_probs[s][i] = -1.0;
            }
        }
    }
    for (int s = 0; s < nSamples; s++) {
        if (NULL != this->base_qScores) {
            this->base_qScores[s] = (int*)realloc(this->base_qScores[s], this->_nBasesPerSample * sizeof(int));
        }
        if (NULL != this->adj_base_qScores) {
            this->adj_base_qScores[s] = (int*)realloc(this->adj_base_qScores[s], this->_nBasesPerSample * sizeof(int));
        }
        this->bases[s] = (int*)realloc(this->bases[s], this->_nBasesPerSample * sizeof(int));
        for (int i = old_nBasesPerSample;i < this->_nBasesPerSample;++i) {
            if (NULL != this->base_qScores) {
                this->base_qScores[s][i] = -1;
            }
            if (NULL != this->adj_base_qScores) {
                this->adj_base_qScores[s][i] = -1;
            }
            this->bases[s][i] = -1;
        }
    }
}

// @return MACROS:CONSTANTS:SIGNALS
// return signals:
// GVCF_NO_WRITE           do not write anything 
//                           i.e. current sim->rec was either: 
//                           1) added to an existing block, or
//                           2) is being held hostage in memory to be the
//                           beginning of a new block
// GVCF_FLUSH_BLOCK        flush the existing gVCF block
//                           N.B. the function should be called again after
//                           flushing to process the current record
// GVCF_WRITE_SIMREC       write sim->rec as a regular record
// @source mostly based on bcftools/gvcf.c
int prepare_gvcf_block(simRecord* sim, gvcfData* gvcfd) {

    int signal = -1;

    do {

        // ------------------------------------------------------- //
        // -> received no simrec
        // i.e. at the end of the program, calling one last time for flushing
        // any existing gVCF block in memory 
        if (NULL == sim->rec) {

            if (0 == gvcfd->current_block_dpr) {
                // --> no existing block in memory
                // do nothing
                signal = GVCF_NO_WRITE;
                break;
            } else {
                // --> found existing block in memory
                signal = GVCF_FLUSH_BLOCK;
                break;
            }
        }


        // ------------------------------------------------------- //
        // -> received a simrec: sim->rec

        if (0 == gvcfd->current_block_dpr) {
            // --> no existing block in memory

            if (sim->nAllelesObserved != 1) {
                signal = GVCF_WRITE_SIMREC;
                break;
            }

        } else {
            // --> found existing block in memory

            // 1) we should flush the existing gVCF block (before checking dps) if:

            // 1.1) at the end of records (sim->rec == NULL)
            // (already handled before, at the beginning of the function)

            // 1.2) block is broken by a variant site 
            if (sim->nAllelesObserved != 1) {
                signal = GVCF_FLUSH_BLOCK;
                break;
            }

            // 1.3) if the current simrec is on a different chromosome
            // N.B. multi chromosome is not well tested
            if (sim->rec->rid != gvcfd->rid) {
                signal = GVCF_FLUSH_BLOCK;
                break;
            }

            // 1.4) encountered a gap 
            // (only possible if args->explode == 0 and args->gvcf == 1)
            // flush the existing record, then continue with the current record candidate
            // to see if it can be the beginning of a new block
            if (sim->rec->pos > gvcfd->end_pos + 1) {
                signal = GVCF_FLUSH_BLOCK;
                break;
            }

        }

        // ------------------------------------------------------- //
        // ok, simrec seems like a candidate so far
        // now check if it is a candidate based on other criteria as well:

        // init variables
        int32_t min_dp = 0;
        int ret = 0;

        ret = bcf_get_format_int32(sim->hdr, sim->rec, "DP", &gvcfd->tmp, &gvcfd->mtmp);
        ASSERT(ret == sim->nSamples);

        // get the minimum DP value across all samples
        min_dp = gvcfd->tmp[0];
        for (int s = 1;s < sim->nSamples;++s) {
            if (min_dp > gvcfd->tmp[s]) {
                min_dp = gvcfd->tmp[s];
            }
        }

        // find the dp_range that this simrec belongs to
        int32_t dp_range = 0;
        int r = 0;
        for (r = 0; r < gvcfd->_block_dps;++r) {
            if (min_dp < gvcfd->block_dps[r]) {
                break; // for
            }
        }
        dp_range = r;

        if (0 != min_dp) {
            bcf_update_info_int32(sim->hdr, gvcfd->grec, "MIN_DP", &min_dp, 1);
        }



        // ------------------------------------------------------- //
        // -> the simrec dp_range is not valid

        // still at the init value 0:
        // -> dp is too small to have this simrec included in any block
        if (!dp_range) {
            if (0 == gvcfd->current_block_dpr) {
                // --> no existing block in memory
                signal = GVCF_WRITE_SIMREC;
                break;
            } else {
                // --> found existing block in memory
                signal = GVCF_FLUSH_BLOCK;
                break;
            }
        }


        // ------------------------------------------------------- //
        // -> the simrec dp_range is valid 

        if (0 != gvcfd->current_block_dpr) {
            // --> found existing block in memory

            // 1.5) simrec dp_range is different than the dp range of the existing block
            if (gvcfd->current_block_dpr != dp_range) {
                // but simrec dp_range is different than the dp range of the existing block 
                // so we need to flush the existing block in memory first
                // then the loop will run again and we will be able to add the simrec to a new block
                // DRAGON:perf could be improved by keeping this information to avoid looping over dps again after coming back from the flush to process the simrec
                signal = GVCF_FLUSH_BLOCK;
                break;
            }
        }

        // ------------------------------------------------------- //
        // -> simrec is now promoted from being a candidate to being a member 
        // congrats simrec!!

        if (0 == gvcfd->current_block_dpr) {
            // --> no existing block in memory
            // so current simrec is the founder of a new block

            if (gvcfd->dp == NULL) {
                gvcfd->dp = (int32_t*)malloc(sizeof(int32_t) * sim->nSamples);
            }
            memcpy(gvcfd->dp, gvcfd->tmp, sizeof(int32_t) * sim->nSamples);

            gvcfd->npl = bcf_get_format_int32(sim->hdr, sim->rec, "PL", &gvcfd->pl, &gvcfd->mpl);
            gvcfd->nqsum = bcf_get_info_float(sim->hdr, sim->rec, "QS", &gvcfd->qsum, &gvcfd->mqsum);
            gvcfd->ngts = bcf_get_genotypes(sim->hdr, sim->rec, &gvcfd->gts, &gvcfd->mgts);

            ASSERT(sim->fmt_dp_arr != NULL);

            gvcfd->rid = sim->rec->rid;
            gvcfd->start_pos = sim->rec->pos;

            gvcfd->alleles.l = 0;

            kputs(sim->rec->d.allele[0], &gvcfd->alleles);
            for (int i = 1;i < sim->rec->n_allele;++i) {
                kputc(',', &gvcfd->alleles);
                kputs(sim->rec->d.allele[i], &gvcfd->alleles);
            }
            gvcfd->min_dp = min_dp;

            gvcfd->current_block_dpr = dp_range;

        } else {

            // --> we already have a block in memory
            // so current simrec will be a member of the existing block

            DEVASSERT(gvcfd->current_block_dpr == dp_range);

            if (gvcfd->min_dp > min_dp) {
                // if new member's min_dp is smaller than the current block min_dp
                // update it to the new min
                gvcfd->min_dp = min_dp;
            }

            for (int s = 0;s < sim->nSamples;++s) {
                if (gvcfd->dp[s] > sim->fmt_dp_arr[s]) {
                    gvcfd->dp[s] = sim->fmt_dp_arr[s];
                }
            }

            ret = bcf_get_format_int32(sim->hdr, sim->rec, "PL", &gvcfd->tmp, &gvcfd->mtmp);
            if (ret >= 0) {
                if (ret != sim->nSamples * 3) {
                    ERROR("Unexpected number of PL values: %d", ret);
                }

                // summarize the PLs across sites in the block
                // we already know REF=X ALT=<UNOBSERVED>
                // so (3*s + 0) should always be 0
                // use the smallest PL values for REF,ALT and ALT,ALT
                for (int s = 0; s < sim->nSamples; s++) {
                    if (gvcfd->pl[3 * s + 1] > sim->pl_arr[3 * s + 1]) {
                        gvcfd->pl[3 * s + 1] = sim->pl_arr[3 * s + 1];
                        gvcfd->pl[3 * s + 2] = sim->pl_arr[3 * s + 2];
                    } else if ((gvcfd->pl[3 * s + 1] == sim->pl_arr[3 * s + 1]) && (gvcfd->pl[3 * s + 2] > sim->pl_arr[3 * s + 2])) {
                        gvcfd->pl[3 * s + 2] = sim->pl_arr[3 * s + 2];
                    }
                }

            } else {
                gvcfd->npl = 0;
            }

        }

        int32_t* tmp = NULL;
        int32_t ntmp = 0;
        if (1 == bcf_get_info_int32(sim->hdr, sim->rec, "END", &tmp, &ntmp)) {
            NEVER;
            gvcfd->end_pos = tmp[0] - 1; // 1-based to 0-based
        } else {
            gvcfd->end_pos = sim->rec->pos;
        }
        signal = GVCF_NO_WRITE;
        break;


    } while (0);



    if (signal == GVCF_NO_WRITE) {
        return(signal);

    } else if (signal == GVCF_FLUSH_BLOCK) {

        gvcfd->end_pos++; // from 0-based to 1-based position

        bcf_clear1(gvcfd->grec);
        gvcfd->grec->rid = gvcfd->rid;
        gvcfd->grec->pos = gvcfd->start_pos;
        gvcfd->grec->rlen = gvcfd->end_pos - gvcfd->start_pos;
        bcf_update_alleles_str(sim->hdr, gvcfd->grec, gvcfd->alleles.s);

        // if block is bigger than 1bp, add END tag
        if (gvcfd->end_pos - gvcfd->start_pos >= 2) {
            bcf_update_info_int32(sim->hdr, gvcfd->grec, "END", &gvcfd->end_pos, 1);
        } else {
            // TODO check how GATK handles this
            // if block is 1bp, no need to add an END tag
        }

        bcf_update_info_int32(sim->hdr, gvcfd->grec, "MIN_DP", &gvcfd->min_dp, 1);
        if (gvcfd->nqsum > 0) {
            bcf_update_info_float(sim->hdr, gvcfd->grec, "QS", gvcfd->qsum, gvcfd->nqsum);
        }
        if (gvcfd->ngts > 0) {
            bcf_update_genotypes(sim->hdr, gvcfd->grec, gvcfd->gts, gvcfd->ngts);
        }
        if (gvcfd->npl > 0) {
            bcf_update_format_int32(sim->hdr, gvcfd->grec, "PL", gvcfd->pl, gvcfd->npl);
        }
        bcf_update_format_int32(sim->hdr, gvcfd->grec, "DP", gvcfd->dp, sim->nSamples);

        gvcfd->current_block_dpr = 0;
        gvcfd->rid = -1;
        gvcfd->npl = 0;
        gvcfd->nqsum = 0;
        gvcfd->ngts = 0;

        return(signal);

    } else if (signal == GVCF_WRITE_SIMREC) {

        return(signal);

    } else {
        NEVER;
    }


}


gvcfData* gvcfData_init(void) {

    gvcfData* gvcfd = new gvcfData;

    gvcfd->grec = bcf_init();

    gvcfd->alleles = KS_INIT;

    const char* dp_ranges = args->gvcf_dps_str;

    // range arg str reading, copied from bcftools/gvcf.c
    int n = 1;
    const char* ss = dp_ranges;
    while (*ss)
    {
        if (*ss == ',') n++;
        ss++;
    }
    gvcfd->_block_dps = n;
    gvcfd->block_dps = (int*)malloc(sizeof(int) * gvcfd->_block_dps);

    n = 0;
    ss = dp_ranges;
    while (*ss)
    {
        char* se = (char*)ss;
        gvcfd->block_dps[n++] = strtol(ss, &se, 10);
        if (se == ss) {
            return(NULL);
        }
        if (*se == ',' && se[1]) { ss = se + 1; continue; } else if (!*se) break;
        return(NULL);
    }

    for (int i = 0;i < gvcfd->_block_dps;++i) {
        if (gvcfd->block_dps[i] < 1) {
            ERROR("Invalid DP range: %d", gvcfd->block_dps[i]);
        }
    }

    return(gvcfd);
}


void gvcfData_destroy(gvcfData* gvcfd) {
    free(gvcfd->block_dps);
    free(gvcfd->dp);
    free(gvcfd->pl);
    if (gvcfd->tmp != NULL) {
        free(gvcfd->tmp);
    }
    free(gvcfd->qsum);
    free(gvcfd->gts);
    free(gvcfd->alleles.s);
    if (gvcfd->grec != NULL) {
        bcf_destroy(gvcfd->grec);
    }
    delete gvcfd;
}



bcf_tag_t bcf_tags[] = {
    /* *********************************************************************** *
     * INFO		Describe the overall variation at site, one value array
     * 			<Number=X>
     * 				    X	values
     * 				    A	one value per alternate allele
     * 				    R	one value for each possible allele
     * (ref+alt) G	one value for each possible genotype .
     * varies/unknown/unbounded Total size needed: X
     *
     * FORMAT	Describe samples, one value array per sample
     * 			<Number=X>
     * 				    X	values (per sample)
     * 				    A	one value per alternate allele (per
     * sample) R	one value for each possible allele (ref+alt) (per
     * sample) G	one value for each possible genotype (per sample) .
     * varies/unknown/unbounded (per sample) Total size needed: X * nSamples
     *
     * ********************************************************************** */

     /* bcf_tag: [FORMAT/GT]
      */

     [GT] =
         {
             .n = FMT_NUMBER_GT,
             .type = BCF_HT_STR,
             .str = "GT",
             .hdr = "##FORMAT=<ID=GT,Number=1,Type=String,Description="
                    "\"Genotype\">",
         },

    /* bcf_tag: [FORMAT/GL]
     */

    [GL] =
        {
            .n = FMT_NUMBER_G,
            .type = BCF_HT_REAL,
            .str = "GL",
            .hdr = "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype "
                   "likelihood in log10 likelihood ratio format\">",
        },

    /* bcf_tag: [FORMAT/GP]
     */

    [GP] =
        {
            .n = FMT_NUMBER_G,
            .type = BCF_HT_REAL,
            .str = "GP",
            .hdr = "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype "
                   "probabilities\">",
        },

    /* bcf_tag: [FORMAT/PL]
     */

    [PL] =
        {
            .n = FMT_NUMBER_G,
            .type = BCF_HT_INT,
            .str = "PL",
            .hdr = "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-"
                   "scaled genotype likelihoods\">",
        },

    /* bcf_tag: [FORMAT/DP]
     */

    [FMT_DP] =
        {
            .n = FMT_NUMBER_1,
            .type = BCF_HT_INT,
            .str = "DP",
            .hdr = "##FORMAT=<ID=DP,Number=1,Type=Integer,Description="
                   "\"Simulated per-sample read depth\">",
        },

    /* bcf_tag: [INFO/DP]
     */

    [INFO_DP] =
        {
            .n = INFO_NUMBER_1,
            .type = BCF_HT_INT,
            .str = "DP",
            .hdr = "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total "
                   "read depth\">",
        },

    /* bcf_tag: [INFO/QS]
     *
     * Normalized phred-scale quality score sums
     * \def qs_vals[nGenotypes]
     * 		   5  == { A, C, G, T, N (unknown reference) }
     * 		   		 (follows the REF,ALT order at each site)
     *
     * e.g.	site with REF=A	ALT=G
     * 		qs_vals[0] 		= sum of quality scores for As
     * 		qs_vals[1] 		= sum of quality scores for Gs
     * 		qs_vals[2|3|4] 	= -1 (undefined)
     */

    [QS] =
        {
            .n = INFO_NUMBER_R_WITH_NONREF,
            .type = BCF_HT_REAL,
            .str = "QS",
            .hdr = "##INFO=<ID=QS,Number=R,Type=Float,Description=\"Normalized "
                   "phred-score allele quality sum. Auxiliary bcf tag used for "
                   "calling\">",
        },

    /*
     * bcf_tag: [INFO/I16]
     *
     * \def i16_vals[16]
     * 		typically produced by bcftools mpileup
     * 		and necessary for genotype calling using bcftools
     * 			val	description
     * 			1	#reference Q13 bases on the forward strand
     * 			2	#reference Q13 bases on the reverse strand
     * 			3	#non-ref Q13 bases on the forward strand
     * 			4	#non-ref Q13 bases on the reverse strand
     * 			5	sum of reference base qualities
     * 			6	sum of squares of reference base qualities
     * 			7	sum of non-ref base qualities
     * 			8	sum of squares of non-ref base qualities
     * 			9	sum of ref mapping qualities
     * 			10	sum of squares of ref mapping qualities
     * 			11	sum of non-ref mapping qualities
     * 			12	sum of squares of non-ref mapping qualities
     * 			13	sum of tail distance for ref bases
     * 			14	sum of squares of tail distance for ref bases
     * 			15	sum of tail distance for non-ref bases
     * 			16	sum of squares of tail distance for non-ref
     *
     * 	source: https://samtools.sourceforge.net/mpileup.shtml
     *
     * \def i16_vals[16]
     * 		typically produced by bcftools mpileup
     * 		and necessary for genotype calling using bcftools
     * 			val	description
     * 			1	#reference bases on the forward strand
     * 			2	#reference bases on the reverse strand
     * 			3	#non-ref bases on the forward strand
     * 			4	#non-ref bases on the reverse strand
     * 			5	sum of reference base qualities
     * 			6	sum of squares of reference base qualities
     * 			7	sum of non-ref base qualities
     * 			8	sum of squares of non-ref base qualities
     * 			9	sum of ref mapping qualities
     * 			10	sum of squares of ref mapping qualities
     * 			11	sum of non-ref mapping qualities
     * 			12	sum of squares of non-ref mapping qualities
     * 			13	sum of tail distance for ref bases
     * 			14	sum of squares of tail distance for ref bases
     * 			15	sum of tail distance for non-ref bases
     * 			16	sum of squares of tail distance for non-ref
     *
     * 	source: bcftools bam2bcf.h bcf_call_ret2_t->anno
     *
     */

    [I16] =
        {
            .n = INFO_NUMBER_16,
            .type = BCF_HT_REAL,
            .str = "I16",
            .hdr = "##INFO=<ID=I16,Number=16,Type=Float,Description="
                   "\"Auxiliary bcf tag used for calling, see description of "
                   "bcf_callret1_t in bam2bcf.h\">",
        },

    /* bcf_tag: [FORMAT/AD]
     */
    [FMT_AD] =
        {
            .n = FMT_NUMBER_R_WITH_NONREF,
            .type = BCF_HT_INT,
            .str = "AD",
            .hdr = "##FORMAT=<ID=AD,Number=R,Type=Integer,Description="
                   "\"Allelic depths for the REF and ALT alleles\">",
        },

    /* bcf_tag: [FORMAT/ADF]
     */
    [FMT_ADF] =
        {
            .n = FMT_NUMBER_R_WITH_NONREF,
            .type = BCF_HT_INT,
            .str = "ADF",
            .hdr =
                "##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic "
                "depths on the forward strand for the REF and ALT alleles\">",
        },

    /* bcf_tag: [FORMAT/ADR]
     */
    [FMT_ADR] =
        {
            .n = FMT_NUMBER_R_WITH_NONREF,
            .type = BCF_HT_INT,
            .str = "ADR",
            .hdr =
                "##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic "
                "depths on the reverse strand for the REF and ALT alleles\">",
        },

    /* bcf_tag: [INFO/AD]
     */
    [INFO_AD] =
        {
            .n = INFO_NUMBER_R_WITH_NONREF,
            .type = BCF_HT_INT,
            .str = "AD",
            .hdr = "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total "
                   "allelic depths for the REF and ALT alleles\">",
        },

    /* bcf_tag: [INFO/ADF]
     */
    [INFO_ADF] =
        {
            .n = INFO_NUMBER_R_WITH_NONREF,
            .type = BCF_HT_INT,
            .str = "ADF",
            .hdr = "##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Total "
                   "allelic depths on the forward strand for the REF and ALT "
                   "alleles\">",
        },

    /* bcf_tag: [INFO/ADR]
     */
    [INFO_ADR] =
        {
            .n = INFO_NUMBER_R_WITH_NONREF,
            .type = BCF_HT_INT,
            .str = "ADR",
            .hdr = "##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Total "
                   "allelic depths on the reverse strand for the REF and ALT "
                   "alleles\">",
        },

};
