#ifndef __BCF_UTILS__
#define __BCF_UTILS__

#include <htslib/kstring.h> // kstring_t

#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <limits>

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <sys/stat.h>
#include <stdio.h>

#include "shared.h"
#include "version.h"

/* ========================================================================== */
/* /BEGIN/ BCF UTILS ======================================================== */

typedef struct simRecord
{

	bcf_hdr_t *hdr = NULL;

	// used to store bcf record
	// either the record copied from input file or blank record created by -explode 1
	bcf1_t *rec=NULL; 


	// acgt_n_q13_bases	number of bases where quality score >= 13 summed across all individuals
	//
	// \def acgt_n_q13_bases[8]
	// 		acgt_n_q13_bases[A-forward,A-reverse,C-forward,C-reverse,G-forward,G-reverse,T-forward,T-reverse]
	//
	//
	// e.g.
	// 		acgt_n_q13_bases[2] == #C bases that are forward and have quality score >= 13
	int* acgt_n_q13_bases=NULL;

	// acgt_sum_taildist sum of tail lengths
	// \def acgt_sum_taildist[4]
	// 		acgt_sum_taildist[A|C|G|T]
	//
	// e.g.
	// 		acgt_sum_taildist[0] = sum of tail lengths for A bases
	float* acgt_sum_taildist=NULL;
	float* acgt_sum_taildist_sq=NULL;
	float* acgt_sum_qs=NULL;
	float* acgt_sum_qs_sq=NULL;


	kstring_t alleles;

	int nSamples = 0;
	int nSites = 0;


	// Number of alleles 
	// if trimAlts==0; should always bbe equal to 4 
	// else; should be equal to the number of alleles observed at site
	int nAlleles = 0; 

	int nGenotypes = 0; // Number of possible genotypes, ==nGLs ==nGPs ==nPLs

	int nHaplotypes = 0; // nSamples * Ploidy


	// acgt2alleles[i] gives the index of the acgt base in ordered bcf REF+ALT alleles
	// i		index of the base in acgt
	// return	index of the base in ordered alleles
	// e.g.
	// allele			A,C,G,T
	// acgt_index i 	0,1,2,3 <- index in base acgt
	//
	// if site has REF=G, ALT=A,T,C
	//				   0      1,2,3 <- index in alleles
	// acgt2alleles[0(A)] = 1 // base A is at index 1 in alleles
	// acgt2alleles[1(C)] = 3 // base C is at index 3 in alleles
	// acgt2alleles[2(G)] = 0 // base G is at index 0 in alleles
	// acgt2alleles[3(T)] = 2 // base T is at index 2 in alleles
	int* acgt2alleles=NULL;

	int* alleles2acgt=NULL;


	// maximum sizes needed to store enum INFO/FORMAT<Number=bcf_tag_number>
	int* max_size_bcf_tag_number=NULL;

	// current sizes in use for enum INFO/FORMAT<Number=bcf_tag_number>
	int* current_size_bcf_tag_number=NULL;


	// TODO delme
	double *gl_vals = NULL;

	double *mps_depths = NULL;

	int32_t *gt_arr = NULL;

	int32_t *fmt_dp_arr = NULL;
	int32_t *fmt_ad_arr = NULL;
	int32_t *fmt_adf_arr = NULL;
	int32_t *fmt_adr_arr = NULL;
	int32_t *info_ad_arr = NULL;
	int32_t *info_adf_arr = NULL;
	int32_t *info_adr_arr = NULL;

	float *gl_arr = NULL;
	float *gp_arr = NULL;
	int32_t *pl_arr = NULL;

	float *qs_arr = NULL;
	float *i16_arr = NULL;


	simRecord(bcf_hdr_t *in_hdr);
	~simRecord();

	void create_hdr(bcf_hdr_t *in_hdr);

	// get sorted indices of values in acgt_arr (descending)
	// without changing acgt_arr itself
	// set alleles2acgt to new indices
	void set_acgt_alleles_luts(int* acgt_arr);

	void reset_rec_objects();

	void add_tags();

	// sim->set_tag_I16(&acgt_n_q13_bases[0][0], &acgt_sum_taildist[0], &acgt_sum_taildist_sq[0], &acgt_sum_qs[0], &acgt_sum_qs_sq[0]);
	void set_tag_I16();

} simRecord;

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
#define N_ENUM_BCF_TAG_NUMBER 8

enum bcf_tag_number
{
	FMT_NUMBER_1, // 	[0], 1*nSamples
	FMT_NUMBER_GT,//	[1], Ploidy*nSamples, for GT
	FMT_NUMBER_G, // 	[2], nGenotypes*nSamples
	FMT_NUMBER_R, // 	[3], nAlleles*nSamples
	INFO_NUMBER_1, // 	[4], 1
	INFO_NUMBER_G, // 	[5], nGenotypes
	INFO_NUMBER_R, // 	[6], nAlleles
	INFO_NUMBER_16 //	[7], 16
};

// inspired by: bcftools/tag2tag
// @typedef struct bcf_tag_t
// @brief Stores information about a bcf tag
// @field n		bcf tag number (expected number of elements in the tag's array)
//TODO add therest
// @note If you add a new enum bcf_tag (nTags++), make sure to
// -> Update bcf_tag_t information in bcf_tags[new_tag] in bcf_utils.cpp
typedef struct
{
	enum bcf_tag_number n;
	int type; // TODO is this needed?
	const char *str = NULL;
	const char *hdr = NULL;
} bcf_tag_t;

// @enum bcf_tag
// @brief Defines the bcf tags used in this program
// @note If you add a new enum bcf_tag, make sure to
// -> Update bcf_tag_t bcf_tags[] in bcf_utils.cpp
// -> Update extern bcf_tag_t bcf_tags[nElements++] in bcf_utils.h
enum bcf_tag
{
	GT, // 1
	GL, // 2
	GP, // 3
	PL, // 4
	FMT_DP, // 5
	INFO_DP, // 6
	QS, // 7
	I16,// 8
	FMT_AD, // 9 
	FMT_ADF,// 10 
	FMT_ADR, // 11
	INFO_AD, // 12
	INFO_ADF,// 13
	INFO_ADR, // 14
};
extern bcf_tag_t bcf_tags[15];

template <typename T>
T *bcf_tag_alloc(enum bcf_tag t, T init_val, const int size)
{

	ASSERT(size > 0);
	T *arr = (T *)malloc(size * sizeof(T));
	for (int i = 0; i < size; ++i)
	{
		arr[i] = init_val;
	}
	ASSERT(NULL != arr);
	return arr;
}

template <typename T>
T *bcf_tag_alloc_max(enum bcf_tag t, T init_val, simRecord* sim)
{
	int size = sim->max_size_bcf_tag_number[bcf_tags[t].n];
	ASSERT(size > 0);
	T *arr = (T *)malloc(size * sizeof(T));
	for (int i = 0; i < size; ++i)
	{
		arr[i] = init_val;
	}
	ASSERT(NULL != arr);
	return arr;
}

template <typename T>
void bcf_tag_reset(T *arr, enum bcf_tag t, T init_val, simRecord* sim)
{

	const int size = sim->max_size_bcf_tag_number[bcf_tags[t].n];
	ASSERT(size > 0);

	for (int i = 0; i < size; ++i)
	{
		arr[i] = init_val;
	}
	return;
}


/* -> MISC -------------------------------------------------------------------*/

// better alternative: use lut

/// @brief nAlleles2nGenotypes - get the number of possible genotypes assuming ploidy==2
/// @param n	number of alleles
/// @return		number of expected genotypes
/// equivalent to (n * (n+1)) / 2
inline int nAlleles2nGenotypes(const int n)
{
	return (((n * (n + 1)) >> 1));
}

/* ========================================================================== */

#endif // __BCF_UTILS__
