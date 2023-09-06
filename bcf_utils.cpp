#include "bcf_utils.h"

void bcf_tag_set_size(enum bcf_tag t, const int size){
	int bcf_tag_size = bcf_tags[t].n;

	if (1 == bcf_tag_size){
		ERROR("Attempted to set size of a hard-fixed size bcf_tag array.[n=%d]",bcf_tag_size);
	}else if (-1 == bcf_tag_size){ // soft-fixed
		bcf_tags[t].n = size;
		return;
	}else if (-2 == bcf_tag_size){
		ERROR("Attempted to set size of a non-fixed size bcf_tag array.[n=%d]",bcf_tag_size);
	}else if ( bcf_tag_size > 1 ){
		ERROR("Attempted to set size of a hard-fixed size or non-fixed size with updated values bcf_tag array.[n=%d]",bcf_tag_size);
	}else{
		NEVER;
	}
}

bcf_tag_t bcf_tags[] =
{
	/* *********************************************************************** *
	 * INFO		describe the overall variation at site, one value array
	 * INFO=<Number=X>
	 * 				A	one value per alternate allele
	 * 				R	one value for each possible allele (including ref)
	 * 				G	one value for each possible genotype
	 * 				.	varies/unknown/unbounded
	 *
	 * FORMAT	describe samples, one value array per sample
	 * FORMAT=<Number=X>
	 * 				  X values per sample
	 *
	 * ********************************************************************** */

	/* bcf_tag: [FORMAT/DP]
	 */

	[DP]    = { 
		.n    = -1,
		.type = BCF_HT_INT,
		.str  = "DP",
		.hdr  = "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Simulated per-sample read depth\">",
	},

	/* bcf_tag: [FORMAT/GT]
	 */

	[GT]   = { 
		.n    = -1,
		.type = BCF_HT_STR,  
		.str  = "GT", 
		.hdr  = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" ,
	},
	
	/* bcf_tag: [FORMAT/GL]
	 */

	[GL]   = { 
		.n    = -2,
		.type = BCF_HT_REAL, 
		.str  = "GL", 
		.hdr  = "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype likelihood in log10 likelihood ratio format\">",
	},

	/* bcf_tag: [FORMAT/GP]
	 */

	[GP]   = { 
		.n    = -2,
		.type = BCF_HT_REAL, 
		.str  = "GP", 
		.hdr  = "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype probabilities\">" ,
	},

	/* bcf_tag: [FORMAT/PL]
	 */

	[PL]   = { 
		.n    = -2,
		.type = BCF_HT_INT,  
		.str  = "PL", 
		.hdr  = "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">" ,
	},

	/* bcf_tag: [INFO/QS]
	 *
	 * Normalized phred-scale quality score sums
	 * \def qs_vals[5]	
	 * 		   5  == { A, C, G, T, N (unknown reference) }
	 * 		   		 (follows the REF,ALT order at each site)
	 *
	 * e.g.	site with REF=A	ALT=G
	 * 		qs_vals[0] 		= sum of quality scores for As
	 * 		qs_vals[1] 		= sum of quality scores for Gs
	 * 		qs_vals[2|3|4] 	= -1 (undefined)
	 */

	[QS]   = { 
		.n    = 5,
		.type = BCF_HT_REAL,
		.str  = "QS",
		.hdr  = "##INFO=<ID=QS,Number=5,Type=Float,Description=\"Normalized phred-score allele quality sum. Auxiliary bcf_tag used for calling\">",
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
	 *
	 */

	[I16]   = { 
		.n    = 16,
		.type = BCF_HT_REAL, 
		.str  = "I16", 
		.hdr  = "##INFO=<ID=I16,Number=16,Type=Float,Description=\"Auxiliary bcf_tag used for calling, see description of bcf_callret1_t in bam2bcf.h\">",
	},
};

