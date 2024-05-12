#ifndef __BCF_UTILS__
#define __BCF_UTILS__

#include <htslib/kstring.h>  // kstring_t
#include <htslib/vcf.h> // bcf_hdr_t, bcf1_t, etc.

#include "shared.h"

typedef struct gvcfData gvcfData;

/* -> BCF TAGS ---------------------------------------------------------------*/

// @enum bcf_tag_number
// @brief Defines the number of elements in a bcf tag array
// @note If you add a new enum bcf_tag_number (nElements++), make sure to
// -> Update simRecord constructor:
// 		init max_size_bcf_tag_number:
// 			set size to new nElements
//			set init values for new bcf_tag_number
// 		init current_size_bcf_tag_number:
// 			set size to new nElements
//			set init values for new bcf_tag_number
#define N_ENUM_BCF_TAG_NUMBER 10

enum bcf_tag_number {
    FMT_NUMBER_1,               // [0], 1*nSamples, FIXED
    FMT_NUMBER_GT,              // [1], Ploidy*nSamples, for GT, FIXED
    FMT_NUMBER_G,               // [2], nGenotypes*nSamples, VARIABLE_PER_REC
    FMT_NUMBER_R,               // [3], nAllelesObserved*nSamples, VARIABLE_PER_REC
    FMT_NUMBER_R_WITH_NONREF,   // [4], if addNonRef (nAllelesObserved+1)*nSamples, else nAllelesObserved*nSamples, VARIABLE_PER_REC
    INFO_NUMBER_1,              // [5], 1, FIXED
    INFO_NUMBER_G,              // [6], nGenotypes, VARIABLE_PER_REC
    INFO_NUMBER_R,              // [7], nAllelesObserved, VARIABLE_PER_REC
    INFO_NUMBER_R_WITH_NONREF,  // [8], if addNonRef nAllelesObserved+1, else nAllelesObserved, VARIABLE_PER_REC
    INFO_NUMBER_16              // [9], 16, FIXED
};

// inspired by: bcftools/tag2tag
// @typedef struct bcf_tag_t
// @brief Stores information about a bcf tag
// @field n		bcf tag number (expected number of elements in the tag's
// array)
// @field type	bcf tag type (BCF_HT_*), see bcf.h
// @field str		bcf tag string
// @field hdr		bcf tag header
// @note If you add a new enum bcf_tag (nTags++), make sure to
// -> Update bcf_tag_t information in bcf_tags[new_tag] in bcf_utils.cpp
typedef struct {
    enum bcf_tag_number n;
    int type;
    const char* str = NULL;
    const char* hdr = NULL;
} bcf_tag_t;

// @enum bcf_tag
// @brief Defines the bcf tags used in this program
// @note If you add a new enum bcf_tag, make sure to
// -> Update bcf_tag_t bcf_tags[] in bcf_utils.cpp
// -> Update extern bcf_tag_t bcf_tags[nElements++] in bcf_utils.h
enum bcf_tag {
    GT,        // 1
    GL,        // 2
    GP,        // 3
    PL,        // 4
    FMT_DP,    // 5
    INFO_DP,   // 6
    QS,        // 7
    I16,       // 8
    FMT_AD,    // 9
    FMT_ADF,   // 10
    FMT_ADR,   // 11
    INFO_AD,   // 12
    INFO_ADF,  // 13
    INFO_ADR,  // 14
};
extern bcf_tag_t bcf_tags[15];


/* ========================================================================== */

typedef struct simRecord {

    // -> same for all records -------------------------------------------
    // only set once and reused for all records
    int nSamples;
    int nHaplotypes;  // nSamples * SIM_PLOIDY
    int nContigs;


    // -> changes per-record ---------------------------------------------
    int nAlleles;
    int nAllelesObserved;
    int nGenotypes;  // Number of possible genotypes, ==nGLs ==nGPs ==nPLs

    // unobserved allele's index in alleles arr
    int allele_unobserved = -1;

    kstring_t* pileup = NULL;

    kstring_t alleles = KS_INIT;


    // a single bcf record object to be reused at each site
    // allocated at the start of the program
    // for each record to be simulated, either the blank record or the input record
    // is copied into this object
    bcf1_t* rec = NULL;

    gvcfData* gvcfd = NULL;

    bcf_hdr_t* hdr = NULL;
    bcf_hdr_t* truth_hdr = NULL;

    // buffer size for arrays
    int _nBasesPerSample;

    // used iff usePreciseGlError == 1 && preCalc == 0
    // \def base_error_probs[nSamples][nBasesPerSample]
    // if (base_error_probs==NULL); then use args->error_prob_forGl instead
    double** base_error_probs = NULL;

    // used iff preCalc == 0
    // \def base_qScores[nSamples][nBasesPerSample]
    // if (base_qScores==NULL); then use args->qScore instead
    int** base_qScores = NULL;
    int** adj_base_qScores = NULL;

    // \def acgtint_bases[nSamples][nBasesPerSample]
    // int representation of base (A=0,C=1,G=2,T=3)
    int** bases = NULL;


    int32_t* gt_arr = NULL;


    // acgt_n_bases_forI16	number of bases where quality score >= 13 summed across
    // all individuals
    //
    // \def acgt_n_bases_forI16[8]
    // 		acgt_n_bases_forI16[A-forward,A-reverse,C-forward,C-reverse,G-forward,G-reverse,T-forward,T-reverse]
    //
    //
    // e.g.
    // 		acgt_n_bases_forI16[2] == #C bases that are forward and have
    // quality score >= 13
    int* acgt_n_bases_forI16 = NULL;

    // acgt_sum_taildist sum of tail lengths
    // \def acgt_sum_taildist[4]
    // 		acgt_sum_taildist[A|C|G|T]
    //
    // e.g.
    // 		acgt_sum_taildist[0] = sum of tail lengths for A bases
    float* acgt_sum_taildist = NULL;
    float* acgt_sum_taildist_sq = NULL;

    int32_t* acgt_info_ad_arr = NULL;
    int32_t* acgt_fmt_ad_arr = NULL;
    int32_t* acgt_fmt_adf_arr = NULL;
    int32_t* acgt_fmt_adr_arr = NULL;

    // per base qscore sums for each sample
    int32_t* acgt_fmt_qsum_arr = NULL;
    // per base qscore squared sums for each sample
    int32_t* acgt_fmt_qsum_sq_arr = NULL;

    // acgt2alleles[i] gives the index of the acgt base in ordered bcf REF+ALT
    // alleles i		index of the base in acgt return	index of
    // the base in ordered alleles e.g. allele			A,C,G,T
    // acgt_index i 	0,1,2,3 <- index in base acgt
    //
    // if site has REF=G, ALT=A,T,C
    //				   0      1,2,3 <- index in alleles
    // acgt2alleles[0(A)] = 1 // base A is at index 1 in alleles
    // acgt2alleles[1(C)] = 3 // base C is at index 3 in alleles
    // acgt2alleles[2(G)] = 0 // base G is at index 0 in alleles
    // acgt2alleles[3(T)] = 2 // base T is at index 2 in alleles
    int* acgt2alleles = NULL;

    int* alleles2acgt = NULL;

    // maximum sizes needed to store enum INFO/FORMAT<Number=bcf_tag_number>
    int* max_size_bcf_tag_number = NULL;

    // current sizes in use for enum INFO/FORMAT<Number=bcf_tag_number>
    int* current_size_bcf_tag_number = NULL;

    // ---
    // always created, only added if addTYPE==1
    int32_t* fmt_dp_arr = NULL;
    int32_t* info_dp_arr = NULL;

    // \def gl_arr[nSamples*nGenotypes]
    // sample_gl_arr = gl_arr[sample_i*nGenotypes]
    // where nGenotype is variable per record
    float* gl_arr = NULL;
    // ---

    int32_t* pl_arr = NULL;
    float* gp_arr = NULL;
    float* qs_arr = NULL;
    float* i16_arr = NULL;

    int32_t* fmt_ad_arr = NULL;
    int32_t* fmt_adf_arr = NULL;
    int32_t* fmt_adr_arr = NULL;

    int32_t* info_ad_arr = NULL;
    int32_t* info_adf_arr = NULL;
    int32_t* info_adr_arr = NULL;

    simRecord(bcf_hdr_t* in_hdr);
    ~simRecord();

    void set_hdr(bcf_hdr_t* in_hdr);

    void prepare_FMT_NUMBER_R_tags(void);
    void prepare_FMT_NUMBER_G_tags(void);


    void expand_arrays(const int new_size);

    void add_tags();


    // reset the reused arrays by setting elements to initial values
    // these were allocated during simRecord construction, and reused in
    // simulating each record
    void reset_rec_objects() {


        // ACGT: A, C, G, T, <*> (non-ref)
        // assumption: MAX_NALLELES == 5
        this->acgt2alleles[0] = -1;
        this->acgt2alleles[1] = -1;
        this->acgt2alleles[2] = -1;
        this->acgt2alleles[3] = -1;
        this->acgt2alleles[4] = -1;

        this->alleles2acgt[0] = -1;
        this->alleles2acgt[1] = -1;
        this->alleles2acgt[2] = -1;
        this->alleles2acgt[3] = -1;
        this->alleles2acgt[4] = -1;

        this->nAlleles = 0;
        this->nAllelesObserved = 0;

        //TODO check why removing this has no effect
        free(this->gt_arr);
        this->gt_arr = NULL;

        // -> reset INFO_ACGT arrays (size=4) 
        // unroll
        this->acgt_info_ad_arr[0] = 0;
        this->acgt_info_ad_arr[1] = 0;
        this->acgt_info_ad_arr[2] = 0;
        this->acgt_info_ad_arr[3] = 0;

        this->acgt_sum_taildist[0] = 0.0;
        this->acgt_sum_taildist[1] = 0.0;
        this->acgt_sum_taildist[2] = 0.0;
        this->acgt_sum_taildist[3] = 0.0;

        this->acgt_sum_taildist_sq[0] = 0.0;
        this->acgt_sum_taildist_sq[1] = 0.0;
        this->acgt_sum_taildist_sq[2] = 0.0;
        this->acgt_sum_taildist_sq[3] = 0.0;

        // -> reset FMT_ACGT arrays (size=4*nSamples)
        for (int i = 0; i < nSamples * 4; ++i) {
            this->acgt_fmt_ad_arr[i] = 0;
            this->acgt_fmt_adf_arr[i] = 0;
            this->acgt_fmt_adr_arr[i] = 0;
            this->acgt_fmt_qsum_arr[i] = 0;
            this->acgt_fmt_qsum_sq_arr[i] = 0;

        }

        // 2 * base + which_strand
        this->acgt_n_bases_forI16[0] = 0;
        this->acgt_n_bases_forI16[1] = 0;
        this->acgt_n_bases_forI16[2] = 0;
        this->acgt_n_bases_forI16[3] = 0;
        this->acgt_n_bases_forI16[4] = 0;
        this->acgt_n_bases_forI16[5] = 0;
        this->acgt_n_bases_forI16[6] = 0;
        this->acgt_n_bases_forI16[7] = 0;

        int size = -1;


        if (this->fmt_dp_arr != NULL) {
            size = this->max_size_bcf_tag_number[bcf_tags[FMT_DP].n];
            for (int i = 0;i < size;++i) {
                this->fmt_dp_arr[i] = 0;
            }
        }

        if (this->info_dp_arr != NULL) {
            size = this->max_size_bcf_tag_number[bcf_tags[INFO_DP].n];
            for (int i = 0;i < size;++i) {
                this->info_dp_arr[i] = 0;
            }
        }

        // GL PL GP share the same size
        size = this->max_size_bcf_tag_number[bcf_tags[GL].n];
        for (int i = 0;i < size;++i) {
            this->gl_arr[i] = -0.0;

            if (NULL != this->gp_arr) {
                this->gp_arr[i] = MINGP;
            }

            if (NULL != this->pl_arr) {
                this->pl_arr[i] = MAXPL;
            }
        }


        if (this->qs_arr != NULL) {
            size = this->max_size_bcf_tag_number[bcf_tags[QS].n];
            for (int i = 0;i < size;++i) {
                this->qs_arr[i] = 0.0;
            }
        }


        if (this->info_ad_arr != NULL) {
            size = this->max_size_bcf_tag_number[bcf_tags[INFO_AD].n];
            for (int i = 0;i < size;++i) {
                this->info_ad_arr[i] = 0;
            }
        }

        if (this->info_adf_arr != NULL) {
            size = this->max_size_bcf_tag_number[bcf_tags[INFO_ADF].n];
            for (int i = 0;i < size;++i) {
                this->info_adf_arr[i] = 0;
            }
        }

        if (this->info_adr_arr != NULL) {
            size = this->max_size_bcf_tag_number[bcf_tags[INFO_ADR].n];
            for (int i = 0;i < size;++i) {
                this->info_adr_arr[i] = 0;
            }
        }

        if (this->fmt_ad_arr != NULL) {
            size = this->max_size_bcf_tag_number[bcf_tags[FMT_AD].n];
            for (int i = 0;i < size;++i) {
                this->fmt_ad_arr[i] = 0;
            }
        }

        if (this->fmt_adf_arr != NULL) {
            size = this->max_size_bcf_tag_number[bcf_tags[FMT_ADF].n];
            for (int i = 0;i < size;++i) {
                this->fmt_adf_arr[i] = 0;
            }
        }

        if (this->fmt_adr_arr != NULL) {
            size = this->max_size_bcf_tag_number[bcf_tags[FMT_ADR].n];
            for (int i = 0;i < size;++i) {
                this->fmt_adr_arr[i] = 0;
            }
        }

        if (this->i16_arr != NULL) {
            // unroll 
            this->i16_arr[0] = 0.0;
            this->i16_arr[1] = 0.0;
            this->i16_arr[2] = 0.0;
            this->i16_arr[3] = 0.0;
            this->i16_arr[4] = 0.0;
            this->i16_arr[5] = 0.0;
            this->i16_arr[6] = 0.0;
            this->i16_arr[7] = 0.0;
            this->i16_arr[8] = 0.0;
            this->i16_arr[9] = 0.0;
            this->i16_arr[10] = 0.0;
            this->i16_arr[11] = 0.0;
            this->i16_arr[12] = 0.0;
            this->i16_arr[13] = 0.0;
            this->i16_arr[14] = 0.0;
            this->i16_arr[15] = 0.0;
        }


        this->alleles.l = 0;
    }

} simRecord;


/* ========================================================================== */
/* -> GVCF -------------------------------------------------------------------*/


// source: mostly based on gvcf_t in bcftools/gvcf.h
struct gvcfData {

    bcf1_t* grec = NULL;

    int* block_dps = NULL;
    int _block_dps = 0;

    // dp range of the current block in memory
    // values:
    // 0      no block exists in memory
    // X      there is a block in memory with a dp range of X
    int current_block_dpr = 0;

    int32_t* dp = NULL;
    int32_t mdp = 0;
    int32_t* pl = NULL;
    int32_t mpl = 0;
    int32_t npl = 0;
    int32_t* tmp = NULL;
    int32_t mtmp = 0;
    int32_t* gts = NULL;
    int32_t ngts = 0;
    int32_t mgts = 0;
    int32_t nqsum = 0;
    int32_t mqsum = 0;
    float* qsum = NULL;
    int32_t rid = -1;
    int32_t start_pos = -1;
    int32_t end_pos = -1;
    int32_t min_dp = 0;
    kstring_t alleles;


};


// @return
// #define GVCF_NO_WRITE 0
// 0            write nothing, still expanding the block
// #define GVCF_FLUSH_BLOCK 1
// 1            flush the current gVCF block (gvcfd->rec)
//              after flushing, this function is called again to see if the current
//              record can be the beginning of a new block:
//              -> if yes, then the second call will return 0. 
//              -> otherwise the second call will return 2.
// #define GVCF_WRITE_SIMREC 2
// 2            no gVCF block to write; write the current regular record (sim->rec)
int prepare_gvcf_block(simRecord* sim, gvcfData* gvcfd);

gvcfData* gvcfData_init(void);

void gvcfData_destroy(gvcfData* gvcfd);


#endif  // __BCF_UTILS__

