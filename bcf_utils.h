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


typedef struct sim_rec{

	bcf_hdr_t* hdr=NULL;

	bcf1_t* rec=NULL;
	bcf1_t* blank_rec=NULL;

	kstring_t* alleles_str=NULL;

	int nSamples=0;
	int nSites=0;

	int site_i=-1;

	//TODO delme
	double* gl_vals=NULL;

	double* mps_depths=NULL;

	int32_t* dp_arr=NULL;
	int32_t* gt_arr=NULL;

	float* gl_arr=NULL;
	float* gp_arr=NULL;
	int32_t* pl_arr=NULL;

	float* qs_arr=NULL;
	float* i16_arr=NULL;

	sim_rec(bcf_hdr_t* in_hdr);
	~sim_rec();

	void set_hdr(bcf_hdr_t* in_hdr);


}sim_rec;

/* -> BCF TAGS ---------------------------------------------------------------*/

// modified from source: bcftools/tag2tag
typedef struct
{
	int n; // expected number of values in tag
		   // default value:
		   // -1	soft-fixed to nSamples, fixed size for all sites
		   // -2	not fixed, value can vary from site to site based on 
		   // 			site specific values 
		   // >0	hard-fixed with a predefined value that does not require
		   // 			information about the data
	int type;
	const char *str=NULL;
	const char *hdr=NULL;
}bcf_tag_t;

enum bcf_tag{DP, GT, GL, GP, PL, QS, I16};
extern bcf_tag_t bcf_tags[7];


template <typename T> T* bcf_tag_alloc(enum bcf_tag t){
	const int bcf_tag_size = bcf_tags[t].n;
	ASSERT(bcf_tag_size>0);
	T* arr = (T*) malloc(bcf_tag_size * sizeof(T));
	ASSERT(NULL!=arr);
	return arr;
}

template <typename T> T* bcf_tag_alloc(enum bcf_tag t, T init_val){
	const int bcf_tag_size = bcf_tags[t].n;
	ASSERT(bcf_tag_size>0);
	T* arr = (T*) malloc(bcf_tag_size * sizeof(T));
	for(int i=0;i<bcf_tag_size;++i){
		arr[i]=init_val;
	}
	ASSERT(NULL!=arr);
	return arr;
}

template <typename T> void bcf_tag_reset(T* arr, enum bcf_tag t, T init_val, const int size){

	ASSERT(size>0);
	bcf_tags[t].n=size;

	ASSERT(NULL!=arr);
	for(int i=0;i<size;++i){
		arr[i]=init_val;
	}
	return;
}

template <typename T> void bcf_tag_reset(T* arr, enum bcf_tag t, T init_val){

	const int size=bcf_tags[t].n;

	ASSERT(NULL!=arr);
	for(int i=0;i<size;++i){
		arr[i]=init_val;
	}
	return;
}

template <typename T> T* bcf_tag_alloc(enum bcf_tag t, T init_val, const int size){

	int bcf_tag_size = bcf_tags[t].n;

	if(bcf_tag_size!=size){
		bcf_tags[t].n=size;
	}

	ASSERT(size>0);
	T* arr = (T*) malloc(size * sizeof(T));
	for(int i=0;i<size;++i){
		arr[i]=init_val;
	}
	ASSERT(NULL!=arr);
	return arr;
}


void bcf_tag_set_size(enum bcf_tag t, const int size);

/* -> MISC -------------------------------------------------------------------*/


/// @brief nAlleles2nGenotypes - get the number of possible genotypes assuming ploidy==2
/// @param n	number of alleles
/// @return		number of expected genotypes				
/// equivalent to (n * (n+1)) / 2
inline int nAlleles2nGenotypes(const int n){
	return ( ((n * (n+1)) >> 1) );
}




/* ========================================================================== */


#endif // __BCF_UTILS__
