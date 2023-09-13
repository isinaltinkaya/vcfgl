#include "bcf_utils.h"

#include "io.h"

sim_rec::sim_rec(bcf_hdr_t *in_hdr)
{

	this->nSamples = bcf_hdr_nsamples(in_hdr);

	this->set_hdr(in_hdr);

	if (args->mps_depths_fn != NULL)
	{
		this->mps_depths = read_depthsFile(args->mps_depths_fn, this->nSamples);

		fprintf(stderr, "\n");
		for (int sample_i = 0; sample_i < this->nSamples; sample_i++)
		{
			fprintf(stderr, "Individual %d mean per-site depth is set to %f\n", sample_i, this->mps_depths[sample_i]);
		}
	}

	this->gt_arr = NULL; // allocated in-place

	this->dp_arr = bcf_tag_alloc<int32_t>(DP, 0);

	const int MAX_TAGARR_SIZE = this->nSamples * MAX_NGTS;

	this->gl_vals = (double *)malloc(MAX_TAGARR_SIZE * sizeof(double));
	ASSERT(NULL != this->gl_vals);
	for (int i = 0; i < MAX_TAGARR_SIZE; ++i)
	{
		this->gl_vals[i] = -0.0;
	}

	this->gl_arr = bcf_tag_alloc<float>(GL, -0.0, MAX_TAGARR_SIZE);

	if (1 == args->addGP)
	{
		this->gp_arr = bcf_tag_alloc<float>(GL, -0.0, MAX_TAGARR_SIZE);
	}

	if (1 == args->addPL)
	{
		this->pl_arr = bcf_tag_alloc<int32_t>(PL, 0, MAX_TAGARR_SIZE);
	}

	if (1 == args->addQS)
	{
		this->qs_arr = bcf_tag_alloc<float>(QS, 0.0, MAX_NQS);
	}

	if (1 == args->addI16)
	{
		this->i16_arr = bcf_tag_alloc<float>(I16, 0.0);
	}
}

sim_rec::~sim_rec()
{

	bcf_hdr_destroy(this->hdr);

	free(this->gt_arr);
	this->gt_arr = NULL;

	// TODO delme
	free(this->gl_vals);
	this->gl_vals = NULL;

	free(this->gl_arr);
	this->gl_arr = NULL;

	if (1 == args->addGP)
	{
		free(this->gp_arr);
		this->gp_arr = NULL;
	}

	if (1 == args->addPL)
	{
		free(this->pl_arr);
		this->pl_arr = NULL;
	}

	if (1 == args->addQS)
	{
		free(this->qs_arr);
		this->qs_arr = NULL;
	}

	if (1 == args->addI16)
	{
		free(this->i16_arr);
		this->i16_arr = NULL;
	}

	free(this->dp_arr);
	this->dp_arr = NULL;
}

void sim_rec::set_hdr(bcf_hdr_t *in_hdr)
{

	this->hdr = bcf_hdr_dup(in_hdr);

	char *DATE_TAG = NULL;

	ASSERT(asprintf(&DATE_TAG, "##fileDate=%s", args->datetime) > 0);
	ASSERT(0 == bcf_hdr_append(this->hdr, DATE_TAG));
	free(DATE_TAG);
	DATE_TAG = NULL;

	char *SOURCE_TAG = NULL;
	ASSERT(asprintf(&SOURCE_TAG, "##source=%s", args->command) > 0);
	ASSERT(0 == bcf_hdr_append(this->hdr, SOURCE_TAG));
	free(SOURCE_TAG);
	SOURCE_TAG = NULL;

	char *SOURCE_VERSION_TAG;
	ASSERT(asprintf(&SOURCE_VERSION_TAG, "##source=vcfgl version: %s", VCFGL_VERSION) > 0);
	ASSERT(0 == bcf_hdr_append(this->hdr, SOURCE_VERSION_TAG));
	free(SOURCE_VERSION_TAG);
	SOURCE_VERSION_TAG = NULL;

	ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[GL].hdr));

	if (1 == args->addGP)
	{
		ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[GP].hdr));
	}

	if (1 == args->addPL)
	{
		ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[PL].hdr));
	}

	if (-999 != args->mps_depth)
	{
		bcf_tag_set_size(DP, nSamples);
		ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[DP].hdr));

		if (1 == args->addI16)
		{
			ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[I16].hdr));
		}

		if (1 == args->addQS)
		{
			ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[QS].hdr));
		}
	}
}

void bcf_tag_set_size(enum bcf_tag t, const int size)
{
	int bcf_tag_size = bcf_tags[t].n;

	if (1 == bcf_tag_size)
	{
		ERROR("Attempted to set size of a hard-fixed size bcf_tag array.[n=%d]", bcf_tag_size);
	}
	else if (-1 == bcf_tag_size)
	{ // soft-fixed
		bcf_tags[t].n = size;
		return;
	}
	else if (-2 == bcf_tag_size)
	{
		ERROR("Attempted to set size of a non-fixed size bcf_tag array.[n=%d]", bcf_tag_size);
	}
	else if (bcf_tag_size > 1)
	{
		ERROR("Attempted to set size of a hard-fixed size or non-fixed size with updated values bcf_tag array.[n=%d]", bcf_tag_size);
	}
	else
	{
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

		[DP] = {
			.n = -1,
			.type = BCF_HT_INT,
			.str = "DP",
			.hdr = "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Simulated per-sample read depth\">",
		},

		/* bcf_tag: [FORMAT/GT]
		 */

		[GT] = {
			.n = -1,
			.type = BCF_HT_STR,
			.str = "GT",
			.hdr = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
		},

		/* bcf_tag: [FORMAT/GL]
		 */

		[GL] = {
			.n = -2,
			.type = BCF_HT_REAL,
			.str = "GL",
			.hdr = "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype likelihood in log10 likelihood ratio format\">",
		},

		/* bcf_tag: [FORMAT/GP]
		 */

		[GP] = {
			.n = -2,
			.type = BCF_HT_REAL,
			.str = "GP",
			.hdr = "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype probabilities\">",
		},

		/* bcf_tag: [FORMAT/PL]
		 */

		[PL] = {
			.n = -2,
			.type = BCF_HT_INT,
			.str = "PL",
			.hdr = "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">",
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

		[QS] = {
			.n = -2,
			.type = BCF_HT_REAL,
			.str = "QS",
			.hdr = "##INFO=<ID=QS,Number=5,Type=Float,Description=\"Normalized phred-score allele quality sum. Auxiliary bcf tag used for calling\">",
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

		[I16] = {
			.n = 16,
			.type = BCF_HT_REAL,
			.str = "I16",
			.hdr = "##INFO=<ID=I16,Number=16,Type=Float,Description=\"Auxiliary bcf tag used for calling, see description of bcf_callret1_t in bam2bcf.h\">",
		},
};
