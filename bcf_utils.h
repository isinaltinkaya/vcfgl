#ifndef __BCF_UTILS__
#define __BCF_UTILS__

#include "shared.h"


/* ========================================================================== */
/* /BEGIN/ BCF UTILS ======================================================== */

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
}
bcf_tag_t;

enum bcf_tag{DP, GT, GL, GP, PL, QS, I16};
extern bcf_tag_t bcf_tags[7];


template <typename T> T* bcf_tag_alloc(enum bcf_tag t){
	const int bcf_tag_size = bcf_tags[t].n;
	ASSERT(bcf_tag_size>0);
	T* obj = (T*) malloc(bcf_tag_size * sizeof(T));
	ASSERT(NULL!=obj);
	return obj;
}

template <typename T> T* bcf_tag_alloc(enum bcf_tag t, T init_val){
	const int bcf_tag_size = bcf_tags[t].n;
	ASSERT(bcf_tag_size>0);
	T* obj = (T*) malloc(bcf_tag_size * sizeof(T));
	for(int i=0;i<bcf_tag_size;++i){
		obj[i]=init_val;
	}
	ASSERT(NULL!=obj);
	return obj;
}

template <typename T> T* bcf_tag_alloc(enum bcf_tag t, T init_val, const int size){

	int bcf_tag_size = bcf_tags[t].n;

	if(bcf_tag_size!=size){
		bcf_tags[t].n=size;
	}

	ASSERT(size>0);
	T* obj = (T*) malloc(size * sizeof(T));
	for(int i=0;i<size;++i){
		obj[i]=init_val;
	}
	ASSERT(NULL!=obj);
	return obj;
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
