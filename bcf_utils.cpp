#include "bcf_utils.h"

#include "io.h"

simRecord::simRecord(bcf_hdr_t *in_hdr)
{
	
	this->rec = bcf_init();
	this->alleles = KS_INITIALIZE;


	this->nSamples = bcf_hdr_nsamples(in_hdr);


	this->alleles2acgt = (int*) malloc(4 * sizeof(int));
	for(int i=0;i<4;++i){
		this->alleles2acgt[i]=-1;
	}
	this->acgt2alleles = (int*) malloc(4 * sizeof(int));
	for(int i=0;i<4;++i){
		this->acgt2alleles[i]=-1;
	}

	this->acgt_n_q13_bases = (int*) calloc(8, sizeof(int));

	this->acgt_sum_taildist = (float*) malloc(4 * sizeof(float));
	this->acgt_sum_taildist[0] = 0.0;
	this->acgt_sum_taildist[1] = 0.0;
	this->acgt_sum_taildist[2] = 0.0;
	this->acgt_sum_taildist[3] = 0.0;

	this->acgt_sum_taildist_sq = (float*) malloc(4 * sizeof(float));
	this->acgt_sum_taildist_sq[0] = 0.0;
	this->acgt_sum_taildist_sq[1] = 0.0;
	this->acgt_sum_taildist_sq[2] = 0.0;
	this->acgt_sum_taildist_sq[3] = 0.0;

	this->acgt_sum_qs = (float*) malloc(4 * sizeof(float));
	this->acgt_sum_qs[0] = 0.0;
	this->acgt_sum_qs[1] = 0.0;
	this->acgt_sum_qs[2] = 0.0;
	this->acgt_sum_qs[3] = 0.0;

	this->acgt_sum_qs_sq = (float*) malloc(4 * sizeof(float));
	this->acgt_sum_qs_sq[0] = 0.0;
	this->acgt_sum_qs_sq[1] = 0.0;
	this->acgt_sum_qs_sq[2] = 0.0;
	this->acgt_sum_qs_sq[3] = 0.0;


	fprintf(stderr, "Number of samples: %i\n", this->nSamples);
	fprintf(stderr, "Number of contigs: %d\n", in_hdr->n[BCF_DT_CTG]);


	// set the maximum size of each bcf tag array
	// then init all current tag array sizes to max
	// if trimming not used then these are not modified so they remain equal to max
	this->max_size_bcf_tag_number = (int*) malloc(N_ENUM_BCF_TAG_NUMBER*sizeof(int));
	this->current_size_bcf_tag_number = (int*) malloc(N_ENUM_BCF_TAG_NUMBER*sizeof(int));

	this->nHaplotypes=this->nSamples * SIM_PLOIDY;
	this->max_size_bcf_tag_number[FMT_NUMBER_GT]=this->nHaplotypes;
	this->current_size_bcf_tag_number[FMT_NUMBER_GT]=this->max_size_bcf_tag_number[FMT_NUMBER_GT];

	this->max_size_bcf_tag_number[FMT_NUMBER_1]=this->nSamples;
	// fixed size for all sites
	this->current_size_bcf_tag_number[FMT_NUMBER_1]=this->max_size_bcf_tag_number[FMT_NUMBER_1]; 

	this->max_size_bcf_tag_number[FMT_NUMBER_G]=this->nSamples * MAX_NGTS;
	// variable size depending on the number of possible genotypes 
	this->current_size_bcf_tag_number[FMT_NUMBER_G]=this->max_size_bcf_tag_number[FMT_NUMBER_G]; 

	this->max_size_bcf_tag_number[FMT_NUMBER_R]=MAX_NALLELES*this->nSamples;
	// variable size depending on the number of possible alleles at each site 
	this->current_size_bcf_tag_number[FMT_NUMBER_R]=this->max_size_bcf_tag_number[FMT_NUMBER_R];

	this->max_size_bcf_tag_number[INFO_NUMBER_1]=1;
	// fixed size for all sites
	this->current_size_bcf_tag_number[INFO_NUMBER_1]=1;

	this->max_size_bcf_tag_number[INFO_NUMBER_G]=this->nSamples * MAX_NGTS;
	this->current_size_bcf_tag_number[INFO_NUMBER_G]=this->max_size_bcf_tag_number[INFO_NUMBER_G];

	this->max_size_bcf_tag_number[INFO_NUMBER_R]=MAX_NALLELES;
	// variable size depending on the number of possible alleles at each site
	this->current_size_bcf_tag_number[INFO_NUMBER_R]=this->max_size_bcf_tag_number[INFO_NUMBER_R]; 

	this->max_size_bcf_tag_number[INFO_NUMBER_16]=16;
	// fixed size for all sites
	this->current_size_bcf_tag_number[INFO_NUMBER_16]=16; 

	this->create_hdr(in_hdr);

	if (args->mps_depths_fn != NULL)
	{
		this->mps_depths = read_depthsFile(args->mps_depths_fn, this->nSamples);

		fprintf(stderr, "\n");
		if (args->verbose>0){
			for (int sample_i = 0; sample_i < this->nSamples; sample_i++)
			{
				fprintf(stderr, "Individual \'%s\' (index:%d) mean per-site depth is set to %f\n", bcf_hdr_int2id(this->hdr, BCF_DT_SAMPLE, sample_i), sample_i, this->mps_depths[sample_i]);
			}
		}
	}

	this->gt_arr = NULL; // allocated in-place



	this->gl_vals = (double *)malloc(this->max_size_bcf_tag_number[bcf_tags[GL].n] * sizeof(double));
	ASSERT(NULL != this->gl_vals);
	for (int i = 0; i < this->max_size_bcf_tag_number[bcf_tags[GL].n]; i++)
	{
		this->gl_vals[i] = -0.0;
	}


	// always created, only added if addTYPE==1
	this->gl_arr = bcf_tag_alloc_max<float>(GL, -0.0, this);
	this->fmt_dp_arr = bcf_tag_alloc_max<int32_t>(FMT_DP, 0, this);
	this->info_ad_arr = bcf_tag_alloc_max<int32_t>(INFO_AD, 0, this);

	if (1 == args->addGP)
	{
		this->gp_arr = bcf_tag_alloc_max<float>(GP, -0.0, this);
	}

	if (1 == args->addPL)
	{
		this->pl_arr = bcf_tag_alloc_max<int32_t>(PL, 0, this);
	}

	if (1 == args->addQS)
	{
		this->qs_arr = bcf_tag_alloc_max<float>(QS, 0.0, this);
	}

	if (1 == args->addI16)
	{
		this->i16_arr = bcf_tag_alloc_max<float>(I16, 0.0, this);
	}

	if (1 == args->addFormatAD){
		this->fmt_ad_arr = bcf_tag_alloc_max<int32_t>(FMT_AD, 0, this);
	}

	if (1 == args->addFormatADF){
		this->fmt_adf_arr = bcf_tag_alloc_max<int32_t>(FMT_ADF, 0, this);
	}

	if (1 == args->addFormatADR){
		this->fmt_adr_arr = bcf_tag_alloc_max<int32_t>(FMT_ADR, 0, this);
	}


}

simRecord::~simRecord()
{

	bcf_destroy(this->rec);

	free(this->alleles.s);
	this->alleles.s = NULL;
	this->alleles.m = 0;
	this->alleles.l = 0;

	free(this->acgt_n_q13_bases);
	this->acgt_n_q13_bases = NULL;

	free(this->acgt_sum_taildist);
	this->acgt_sum_taildist = NULL;

	free(this->acgt_sum_taildist_sq);
	this->acgt_sum_taildist_sq = NULL;

	free(this->acgt_sum_qs);
	this->acgt_sum_qs = NULL;

	free(this->acgt_sum_qs_sq);
	this->acgt_sum_qs_sq = NULL;


	bcf_hdr_destroy(this->hdr);

	free(this->gt_arr);
	this->gt_arr = NULL;

	// TODO delme
	free(this->gl_vals);
	this->gl_vals = NULL;

	free(this->gl_arr);
	this->gl_arr = NULL;

	free(this->fmt_dp_arr);
	this->fmt_dp_arr = NULL;

	free(this->info_ad_arr);
	this->info_ad_arr = NULL;


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

	if (1 == args->addFormatAD){
		free(this->fmt_ad_arr);
		this->fmt_ad_arr = NULL;
	}

	if (1 == args->addFormatADF){
		free(this->fmt_adf_arr);
		this->fmt_adf_arr = NULL;
	}

	if (1 == args->addFormatADR){
		free(this->fmt_adr_arr);
		this->fmt_adr_arr = NULL;
	}


	free(this->max_size_bcf_tag_number);
	this->max_size_bcf_tag_number = NULL;

	free(this->current_size_bcf_tag_number);
	this->current_size_bcf_tag_number = NULL;

	if (NULL != this->mps_depths)
	{
		free(this->mps_depths);
		this->mps_depths = NULL;
	}

	free(this->alleles2acgt);
	this->alleles2acgt=NULL;
	free(this->acgt2alleles);
	this->acgt2alleles=NULL;



}

void simRecord::add_tags(){

	if (1 == args->addInfoAD)
	{
		ASSERT(0 == (bcf_update_info_int32(this->hdr,this->rec, "AD", this->info_ad_arr, this->current_size_bcf_tag_number[bcf_tags[INFO_AD].n])));
	}

	if (1 == args->addFormatDP)
	{
		ASSERT(0 == (bcf_update_format_int32(this->hdr,this->rec, "DP", this->fmt_dp_arr, this->current_size_bcf_tag_number[bcf_tags[FMT_DP].n])));
	}

	if (1 == args->addFormatAD)
	{
		ASSERT(0 == (bcf_update_format_int32(this->hdr,this->rec, "AD", this->fmt_ad_arr, this->current_size_bcf_tag_number[bcf_tags[FMT_AD].n])));
	}
	if (1 == args->addFormatADF)
	{
		ASSERT(0 == (bcf_update_format_int32(this->hdr,this->rec, "ADF", this->fmt_adf_arr, this->current_size_bcf_tag_number[bcf_tags[FMT_ADF].n])));
	}
	if (1 == args->addFormatADR)
	{
		ASSERT(0 == (bcf_update_format_int32(this->hdr,this->rec, "ADR", this->fmt_adr_arr,this->current_size_bcf_tag_number[bcf_tags[FMT_ADR].n])));
	}


	if (1 == args->addGL)
	{
		ASSERT(0 == (bcf_update_format_float(this->hdr,this->rec, "GL", this->gl_arr, this->current_size_bcf_tag_number[bcf_tags[GL].n])));
	}

	if (1 == args->addPL)
	{
		ASSERT(0 == (bcf_update_format_int32(this->hdr, this->rec, "PL", this->pl_arr, this->current_size_bcf_tag_number[bcf_tags[PL].n])));
	}

	if (1 == args->addGP)
	{
		ASSERT(0 == (bcf_update_format_float(this->hdr, this->rec, "GP", this->gp_arr, this->current_size_bcf_tag_number[bcf_tags[GP].n])));
	}

	if (1 == args->addI16)
	{

		ASSERT(0 == (bcf_update_info_float(this->hdr, this->rec, "I16", this->i16_arr, 16)));
	}

	// [BEGIN ADD QS TAG] ------------------------------------------------------- //
	if (1 == args->addQS)
	{
		// already reordered while setting sim->qs_arr
		ASSERT(0 == (bcf_update_info_float(this->hdr, this->rec, "QS", this->qs_arr, this->current_size_bcf_tag_number[bcf_tags[QS].n])));
	}



}

void simRecord::reset_rec_objects(){

	// reset the reused arrays by setting elements to initial values
	// these were allocated during simRecord construction, and reused in simulating each record

	this->nAlleles=0;

	this->acgt_n_q13_bases[0]=0;
	this->acgt_n_q13_bases[1]=0;
	this->acgt_n_q13_bases[2]=0;
	this->acgt_n_q13_bases[3]=0;
	this->acgt_n_q13_bases[4]=0;
	this->acgt_n_q13_bases[5]=0;
	this->acgt_n_q13_bases[6]=0;
	this->acgt_n_q13_bases[7]=0;

	this->acgt_sum_taildist[0] = 0.0;
	this->acgt_sum_taildist[1] = 0.0;
	this->acgt_sum_taildist[2] = 0.0;
	this->acgt_sum_taildist[3] = 0.0;

	this->acgt_sum_taildist_sq[0] = 0.0;
	this->acgt_sum_taildist_sq[1] = 0.0;
	this->acgt_sum_taildist_sq[2] = 0.0;
	this->acgt_sum_taildist_sq[3] = 0.0;

	this->acgt_sum_qs[0] = 0.0;
	this->acgt_sum_qs[1] = 0.0;
	this->acgt_sum_qs[2] = 0.0;
	this->acgt_sum_qs[3] = 0.0;

	this->acgt_sum_qs_sq[0] = 0.0;
	this->acgt_sum_qs_sq[1] = 0.0;
	this->acgt_sum_qs_sq[2] = 0.0;
	this->acgt_sum_qs_sq[3] = 0.0;


	// always created, only added if addTYPE==1
	bcf_tag_reset<int32_t>(this->fmt_dp_arr, FMT_DP, 0, this);
	bcf_tag_reset<float>(this->gl_arr, GL, MINGL, this);
	bcf_tag_reset<int32_t>(this->info_ad_arr, INFO_AD, 0, this);

	if (1 == args->addFormatAD)
	{
		bcf_tag_reset<int32_t>(this->fmt_ad_arr, FMT_AD, 0, this);
	}
	if (1 == args->addFormatADF)
	{
		bcf_tag_reset<int32_t>(this->fmt_adf_arr, FMT_ADF, 0, this);
	}
	if (1 == args->addFormatADR)
	{
		bcf_tag_reset<int32_t>(this->fmt_adr_arr, FMT_ADR, 0, this);
	}

	if (1 == args->addQS)
	{
		bcf_tag_reset<float>(this->qs_arr, QS, 0, this);
	}

	if (1 == args->addGP)
	{
		bcf_tag_reset<float>(this->gp_arr, GP, MINGP, this);
	}

	if (1 == args->addPL)
	{
		bcf_tag_reset<int32_t>(this->pl_arr, PL, MINPL , this);
	}

	this->alleles.l = 0;

}

void simRecord::create_hdr(bcf_hdr_t *in_hdr)
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

	if (1 == args->addGL){
		ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[GL].hdr));
	}

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
		if(1==args->addFormatDP){
			ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[FMT_DP].hdr));
		}

		if (1 == args->addI16)
		{
			ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[I16].hdr));
		}

		if (1 == args->addQS)
		{
			ASSERT(0 == bcf_hdr_append(this->hdr, bcf_tags[QS].hdr));
		}
	}


	if (1 == args->addFormatAD){
		ASSERT(0==bcf_hdr_append(this->hdr, bcf_tags[FMT_AD].hdr));
	}
	if (1 == args->addFormatADF){
		ASSERT(0==bcf_hdr_append(this->hdr, bcf_tags[FMT_ADF].hdr));
	}
	if (1 == args->addFormatADR){
		ASSERT(0==bcf_hdr_append(this->hdr, bcf_tags[FMT_ADR].hdr));
	}

	if (1 == args->addInfoAD){
		ASSERT(0==bcf_hdr_append(this->hdr, bcf_tags[INFO_AD].hdr));
	}


}


void simRecord::set_acgt_alleles_luts(int* acgt_arr)
{
	// DEVPRINT("%c:%d %c:%d %c:%d %c:%d","ACGT"[0],acgt_arr[0],"ACGT"[1],acgt_arr[1],"ACGT"[2],acgt_arr[2],"ACGT"[3],acgt_arr[3]);
	ASSERT(NULL!=acgt_arr);
	ASSERT(NULL!=alleles2acgt);
	ASSERT(acgt_arr[0]>-1);
	ASSERT(acgt_arr[1]>-1);
	ASSERT(acgt_arr[2]>-1);
	ASSERT(acgt_arr[3]>-1);

	int* ptr[5];
	// TODO make this ptr[4] then check if tests still pass
	int* tmp;
	
	for (int i = 0; i < 4; i++)
	{
		ptr[i] = &acgt_arr[i];
	}
	// insertion sort
	for (int i = 1; i < 4; i++)
	{
		for (int j = i; j > 0 && *ptr[j] < *ptr[j - 1]; j--)
		{
			tmp = ptr[j], ptr[j] = ptr[j - 1], ptr[j - 1] = tmp;
		}
	}
	for(int i=0;i<4;++i){
		// indices lut descending order
		this->alleles2acgt[i]=ptr[3-i]-acgt_arr;
	}
	for(int i=0;i<4;++i){
		this->acgt2alleles[this->alleles2acgt[i]]=i;
	}
	// DEVPRINT("-> ordered %c:%d %c:%d %c:%d %c:%d","ACGT"[alleles2acgt[0]],acgt_arr[alleles2acgt[0]],"ACGT"[alleles2acgt[1]],acgt_arr[alleles2acgt[1]],"ACGT"[alleles2acgt[2]],acgt_arr[alleles2acgt[2]],"ACGT"[alleles2acgt[3]],acgt_arr[alleles2acgt[3]]);
	// DEVPRINT("-> acgt %c:%d %c:%d %c:%d %c:%d","ACGT"[0],acgt_arr[0],"ACGT"[1],acgt_arr[1],"ACGT"[2],acgt_arr[2],"ACGT"[3],acgt_arr[3]);
	


	return;
}

void simRecord::set_tag_I16(void){

	if(1==args->addI16){

	// ------------------------------------------------------- //
	// I16 REF fields
	// ------------------------------------------------------- //
	int refb=this->alleles2acgt[0];

	// 1   #reference Q13 bases on the forward strand
	this->i16_arr[0] = acgt_n_q13_bases[(refb*2)+SIM_FORWARD_STRAND];

	// 2   #reference Q13 bases on the reverse strand
	this->i16_arr[1] = acgt_n_q13_bases[(refb*2)+SIM_REVERSE_STRAND];

	// 5   sum of reference base qualities
	this->i16_arr[4] = acgt_sum_qs[refb];

	// 6   sum of squares of reference base qualities
	this->i16_arr[5] = acgt_sum_qs_sq[refb];

	// 13  sum of tail distance for ref bases
	this->i16_arr[12] = acgt_sum_taildist[refb];

	// 14  sum of squares of tail distance for ref bases
	this->i16_arr[13] = acgt_sum_taildist_sq[refb];

	// ------------------------------------------------------- //
	// I16 NON-REF fields
	// ------------------------------------------------------- //
	// start from 1 to exclude the reference allele
	int b=-1;
	ASSERT(this->nAlleles>1);
	for (int a = 1; a < this->nAlleles; ++a)
	{
		b=this->alleles2acgt[a];

		// 3   #non-ref Q13 bases on the forward strand
		this->i16_arr[2] += acgt_n_q13_bases[(b*2)+SIM_FORWARD_STRAND];

		// 4   #non-ref Q13 bases on the reverse strand
		this->i16_arr[3] += acgt_n_q13_bases[(b*2)+SIM_REVERSE_STRAND];

		// 7   sum of non-ref base qualities
		this->i16_arr[6] += acgt_sum_qs[b];

		// 8   sum of squares of non-ref base qualities
		this->i16_arr[7] += acgt_sum_qs_sq[b];

		// 15  sum of tail distance for non-ref bases
		this->i16_arr[14] += acgt_sum_taildist[b];

		// 16  sum of squares of tail distance for non-ref
		this->i16_arr[15] += acgt_sum_taildist_sq[b];


	}

	// ------------------------------------------------------- //
	// I16 fields that are currently not simulated
	// and are set to 0
	// ------------------------------------------------------- //
	// TODO

	// 9   sum of ref mapping qualities
	this->i16_arr[8] = 0;

	// 10  sum of squares of ref mapping qualities
	this->i16_arr[9] = 0;

	// 11  sum of non-ref mapping qualities
	this->i16_arr[10] = 0;

	// 12  sum of squares of non-ref mapping qualities
	this->i16_arr[11] = 0;

}

}

bcf_tag_t bcf_tags[] =
	{
		/* *********************************************************************** *
		 * INFO		Describe the overall variation at site, one value array
		 * 			<Number=X>
		 * 				    X	values 
		 * 				    A	one value per alternate allele
		 * 				    R	one value for each possible allele (ref+alt)
		 * 				    G	one value for each possible genotype
		 * 				    .	varies/unknown/unbounded
		 * 			Total size needed: X
		 * 
		 * FORMAT	Describe samples, one value array per sample
		 * 			<Number=X>
		 * 				    X	values (per sample)
		 * 				    A	one value per alternate allele (per sample)
		 * 				    R	one value for each possible allele (ref+alt) (per sample)
		 * 				    G	one value for each possible genotype (per sample)
		 * 				    .	varies/unknown/unbounded (per sample)
		 * 			Total size needed: X * nSamples
		 *
		 * ********************************************************************** */


		/* bcf_tag: [FORMAT/GT]
		 */

		[GT] = {
			.n = FMT_NUMBER_GT,
			.type = BCF_HT_STR,
			.str = "GT",
			.hdr = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
		},

		/* bcf_tag: [FORMAT/GL]
		 */

		[GL] = {
			.n = FMT_NUMBER_G,
			.type = BCF_HT_REAL,
			.str = "GL",
			.hdr = "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype likelihood in log10 likelihood ratio format\">",
		},

		/* bcf_tag: [FORMAT/GP]
		 */

		[GP] = {
			.n = FMT_NUMBER_G,
			.type = BCF_HT_REAL,
			.str = "GP",
			.hdr = "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype probabilities\">",
		},

		/* bcf_tag: [FORMAT/PL]
		 */

		[PL] = {
			.n = FMT_NUMBER_G,
			.type = BCF_HT_INT,
			.str = "PL",
			.hdr = "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">",
		},

		/* bcf_tag: [FORMAT/FMT_DP]
		 */

		[FMT_DP] = {
			.n = FMT_NUMBER_1,
			.type = BCF_HT_INT,
			.str = "DP",
			.hdr = "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Simulated per-sample read depth\">",
		},

		[INFO_DP] = {
			.n = INFO_NUMBER_1,
			.type = BCF_HT_INT,
			.str = "DP",
			.hdr = "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">",
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
			.n = INFO_NUMBER_R,
			.type = BCF_HT_REAL,
			.str = "QS",
			.hdr = "##INFO=<ID=QS,Number=R,Type=Float,Description=\"Normalized phred-score allele quality sum. Auxiliary bcf tag used for calling\">",
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
			.n = INFO_NUMBER_16,
			.type = BCF_HT_REAL,
			.str = "I16",
			.hdr = "##INFO=<ID=I16,Number=16,Type=Float,Description=\"Auxiliary bcf tag used for calling, see description of bcf_callret1_t in bam2bcf.h\">",
		},

		/* bcf_tag: [FORMAT/AD]
		 */
		[FMT_AD] = {
			.n = FMT_NUMBER_R,
			.type = BCF_HT_INT,
			.str = "AD",
			.hdr = "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the REF and ALT alleles\">",
		},

		/* bcf_tag: [FORMAT/ADF]
		 */
		[FMT_ADF] = {
			.n = FMT_NUMBER_R,
			.type = BCF_HT_INT,
			.str = "ADF",
			.hdr = "##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on the forward strand for the REF and ALT alleles\">",
		},

		/* bcf_tag: [FORMAT/ADR]
		 */
		[FMT_ADR] = {
			.n = FMT_NUMBER_R,
			.type = BCF_HT_INT,
			.str = "ADR",
			.hdr = "##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on the reverse strand for the REF and ALT alleles\">",
		},

		/* bcf_tag: [INFO/AD]
		 */
		[INFO_AD] = {
			.n = INFO_NUMBER_R,
			.type = BCF_HT_INT,
			.str = "AD",
			.hdr = "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total allelic depths for the REF and ALT alleles\">",
		},

		/* bcf_tag: [INFO/ADF]
		 */
		[INFO_ADF] = {
			.n = INFO_NUMBER_R,
			.type = BCF_HT_INT,
			.str = "ADF",
			.hdr = "##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Total allelic depths on the forward strand for the REF and ALT alleles\">",
		},

		/* bcf_tag: [INFO/ADR]
		 */
		[INFO_ADR] = {
			.n = INFO_NUMBER_R,
			.type = BCF_HT_INT,
			.str = "ADR",
			.hdr = "##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Total allelic depths on the reverse strand for the REF and ALT alleles\">",
		},


};
