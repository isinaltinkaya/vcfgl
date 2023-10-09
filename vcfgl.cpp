/*
 * vcfgl: gl simulator for vcf/bcf files
 *
 * isinaltinkaya
 *
 *
 */


#include "shared.h"
#include "random_generator.h"
#include "estimator.h"
#include "io.h"
#include "bcf_utils.h"

#include <htslib/thread_pool.h>

argStruct *args;

int (*get_strand)(void);

// a negative value is returned iff 1==args->rmInvarSites
// negative return values indicate skip the site
// value of the negative return value indicate the reason for skipping
// -1	observed only 1 allele at site
// -2	no alleles were observed at site
int (*simulate_record)(simRecord *sim);

int set_strand_to_forward(void){
	return(SIM_FORWARD_STRAND);
}

int simulate_record_values(simRecord *sim)
{

	bcf1_t *rec = sim->rec;

	int sOffset = 0; // sample offset

	double r_error_prob = args->error_rate; // init val
	const int nSamples = sim->nSamples;

	int acgt_fmt_ad_arr[nSamples * 4];
	for (int i = 0; i < nSamples * 4; ++i)
	{
		acgt_fmt_ad_arr[i] = 0;
	}

	sim->reset_rec_objects();

	int n_sim_reads = -1;

	int qScore = -1;

	int which_strand, which_haplo = -1;

	int true_base = -1; // true base (no error)
	int r_base = -1;	// simulated base

	int32_t ngt_arr = 0;
	int ngt = bcf_get_genotypes(sim->hdr, rec, &sim->gt_arr, &ngt_arr);
	if (ngt <= 0)
	{
		ERROR("Could not find GT tag.");
	}


	// per base qscore sums for each sample
	// 	where alleles = {A, C, G, T} (thus [4])
	int samples_acgt_qs[4][nSamples];
	for (int j = 0; j < 4; ++j)
	{
		for (int i = 0; i < nSamples; ++i)
		{
			samples_acgt_qs[j][i] = 0.0;
		}
	}

	// A,C,G,T,N
	// int seens[5] = {0, 1, 2, 3, 4};

	int tail_dist = 0;

	for (int sample_i = 0; sample_i < nSamples; sample_i++)
	{

		if (sim->mps_depths != NULL)
		{
			n_sim_reads = Poisson(sim->mps_depths[sample_i]);
		}
		else
		{
			n_sim_reads = Poisson(args->mps_depth);
		}

		if (n_sim_reads == 0)
		{
			sim->fmt_dp_arr[sample_i] = 0;
		}
		else
		{

			int32_t *gt_ptr = sim->gt_arr + sample_i * SIM_PLOIDY;

			double *sample_gls = sim->gl_vals + sample_i * MAX_NGTS;

			int bin_gts[2] = {0};

			// TODO make it possible to use alt ref other than 0 1 

			// binary input genotypes from simulated input
			for (int i = 0; i < SIM_PLOIDY; i++)
			{
				bin_gts[i] = bcf_gt_allele(gt_ptr[i]);
				

				// TODO gets base form
				// use bcf_unpack(rec,BCF_UN_STR)
				// DEVPRINT("%d",*rec->d.allele[bcf_gt_allele(gt_ptr[0])]);
				// if 0 1 (currently supported notation), this should be 0 or 1
				
				// use bit shifting to check if bin_gt[n] is 0 or 1
				if ((bin_gts[i] >> 1) != 0)
				{
					fprintf(stderr, "ERROR:\n\nbin_gts[%d]: Genotype %d not supported.\n", i, bin_gts[i]);
					exit(1);
				}
			}

			// int sim_strands[n_sim_reads] = {-1};

			// for(int i=0;i<10;i++){
				// ASSERT(sample_gls[i]==-0.0);
			// }

			for (int i = 0; i < n_sim_reads; i++)
			{

				(drand48() < 0.5) ? which_haplo = 0 : which_haplo = 1;

				r_base = true_base = bin_gts[which_haplo];

				if (drand48() < args->error_rate)
				{
					while ((r_base = (floor(4 * drand48()))) == true_base);
				}

				if (1 == args->error_qs)
				{
					r_error_prob = args->betaSampler->sample();
				}

				if (0==r_error_prob){
					qScore = CAP_BASEQ;
				}else{
					if (0 == args->error_qs)
					{
						qScore = args->error_rate_q;
					}
					else
					{
						qScore = floor(-10 * log10(r_error_prob));
					}

				}
				if (qScore > CAP_BASEQ)
				{
					qScore = CAP_BASEQ;
				}

				if (1 == args->platform)
				{

					// [R]
					// qToP<-function(q){10^(-q/10)}
					// 
					// > qToP(2)
					// [1] 0.6309573
					// > qToP(12)
					// [1] 0.06309573
					// > qToP(23)
					// [1] 0.005011872
					// > qToP(37)
					// [1] 0.0001995262

					if (qScore <= 2)
					{
						qScore = 2;
						r_error_prob=0.6309573;
					}
					else if (qScore <= 14)
					{
						qScore = 12;
						r_error_prob=0.06309573;
					}
					else if (qScore <= 30)
					{
						qScore = 23;
						r_error_prob=0.005011872;
					}
					else
					{
						qScore = 37;
						r_error_prob=0.0001995262;
					}

				}


				acgt_fmt_ad_arr[(4 * sample_i) + r_base]++;

				which_strand=get_strand();

				// sim_strands[i] = which_strand;

				if (qScore >= 13)
				{
					sim->acgt_n_q13_bases[(2*r_base) + which_strand]++;
				}

				// fprintf(stdout,"sampled base %d (%c) at site %ld for individual %d with qscore %d\n",r_base,"ACGTN"[r_base],rec->pos+1,sample_i,qScore);

				// ind,site,read,base,error_prob,qscore
				if(sample_i==0){
					fprintf(stdout,"%d,%ld,%d,%c,%f,%d\n",sample_i,rec->pos+1,i,"ACGTN"[r_base],r_error_prob,qScore);
				}

				sim->acgt_sum_qs[r_base] += qScore;
				sim->acgt_sum_qs_sq[r_base] += qs_to_qs2(qScore);

				tail_dist = sample_tail_distance();
				sim->acgt_sum_taildist[r_base] += tail_dist;
				sim->acgt_sum_taildist_sq[r_base] += (tail_dist * tail_dist);

				samples_acgt_qs[r_base][sample_i] += qScore;

				gl_log10(r_base, r_error_prob, sample_gls);
			}

			sim->fmt_dp_arr[sample_i] = n_sim_reads;
			sim->info_dp_arr[0]+=n_sim_reads;
		}

	} // end sample loop
	

	if(0 == sim->info_dp_arr[0]){
	// no reads were simulated for any of the individuals

		VWARN("No alleles were observed at site %ld.", rec->pos + 1);
		if (1 == args->rmInvarSites)
		{
			return (-2);
		}
		if ( 1 == args->trimAlts){
			if(0==args->useUnknownAllele){
				NEVER;
			}
			// REF=<*> ALT=.
			NEVER; //TODO fornow
		}
		// set all requested fields to missing
		sim->nAlleles=4;
		sim->nGenotypes=10;
		ASSERT(0 == (bcf_update_alleles_str(sim->hdr, rec, "A,C,G,T")));
		sim->add_tags_missing();
		return (0);
	}

	int acgt_info_ad_arr[4] = {0};
	sOffset = 0;
	for (int s = 0; s < nSamples; ++s)
	{
		sOffset = s * 4;
		for (int b = 0; b < 4; ++b)
		{
			acgt_info_ad_arr[b] += acgt_fmt_ad_arr[sOffset + b];
		}
	}

	sim->set_acgt_alleles_luts(acgt_info_ad_arr);

	int nUnobservedBases=0;
	for (int b = 0; b < 4; ++b)
	{
		if ( 0 == acgt_info_ad_arr[b])
		{
			if ( 1 == args->trimAlts){
				int tmpa = sim->acgt2alleles[b];
				sim->acgt2alleles[b] = -1;
				sim->alleles2acgt[tmpa] = -1;
			}
			nUnobservedBases++;
		}
	}

	sim->nAlleles = 0;
	for (int a = 0; a < 4; ++a)
	{
		if (-1 == sim->alleles2acgt[a]){
			continue;
		}

		kputc("ACGT"[sim->alleles2acgt[a]], &sim -> alleles);
		sim->nAlleles++;
		if (a != 3 && -1 != sim->alleles2acgt[a + 1])
		{
			kputc(',', &sim->alleles);
		}
	}

	int nObserved = sim->nAlleles - nUnobservedBases;


	// useUnknownAllele
	// trimAlts
	// rmInvarSites
	// 
	// at a site
	// if nAlleles == 0
	// 		observed nothing at all, all individuals are missing
	// 		if rmInvarSites == 1
	// 			return -1; //skip site
	// 		else
	// 			if trimAlts == 1
	// 				never; trimAlts 1 cannot be used with rmInvarSites 1 
	// 			else
	// 				set all requested tags to missing
	// 				set alleles to REF:A ALT:C,G,T; no need to sort and update gts
	//
	// if nAlleles == 1
	// 		observed only one allele at site
	//
	//			if rmInvarSites == 1
	// 				return -1; // skip site
	//			else (rmInvarSites==0)
	////				if trimAlts == 1
	////					trim out the unobserved alleles
	////					if useUnknownAllele == 1
	////						trimAlts 1 cannot be used with useUnknownAllele 1
	////						NEVER;
	// disable trimalts altogether
	////				else (trimAlts==0)
	//			if useUnknownAllele == 1
	//				set ALT to <*> (unobserved)
	//				set nAlleles=2
	//			else (useUnknownAllele==0)
	//				populate non-ref bases ({A,C,G,T} - {REF}, e.g. {C,G,T} for REF=A)
	//				set nAlleles=4

	if (0 == nObserved)
	{
		NEVER;
	}
	else if (1 == nObserved)
	{
		VWARN("Observed only one allele at site %ld.", rec->pos + 1);
		if (1 == args->rmInvarSites)
		{
			return (-1);
		}
		if ( 1 == args->trimAlts){
			// REF=ALLELE, ALT=<*>
			// NEVER; //TODO fornow
		}
	}

	sim->nGenotypes = nAlleles_to_nGenotypes(sim->nAlleles);

	int has_unobserved = 0;


	// TODO if nGenotypes == 1, only one allele was observed, then set the alt to <*> and set nGenotypes to 2
	if (1 == sim->nGenotypes)
	{

		if(1==args->useUnknownAllele){

			// set ALT to <*> (unobserved)

			kputc(',', &sim->alleles);
			kputs("<*>", &sim->alleles);
			sim->nAlleles++;

			has_unobserved = 1;

			sim->nGenotypes = nAlleles_to_nGenotypes(sim->nAlleles);

		}else{
			NEVER;

		}

	}

	sim->current_size_bcf_tag_number[FMT_NUMBER_G] = sim->nSamples * sim->nGenotypes;
	sim->current_size_bcf_tag_number[FMT_NUMBER_R] = sim->nSamples * sim->nAlleles;
	sim->current_size_bcf_tag_number[INFO_NUMBER_G] = sim->nGenotypes;
	sim->current_size_bcf_tag_number[INFO_NUMBER_R] = sim->nAlleles;

	const int nGenotypes = sim->nGenotypes;

	ASSERT(0 == (bcf_update_alleles_str(sim->hdr, rec, sim->alleles.s)));
	ASSERT(sim->nAlleles == rec->n_allele);


	// set sorted values for FMT_NUMBER_R tags
	// 		FMT_AD
	//		FMT_ADF
	//		FMT_ADR
	// if (args->addFormatAD || args->addFormatADF || args->addFormatADR)//TODO
	if (args->addFormatAD)
	{
		sOffset = 0;
		for (int s = 0; s < nSamples; ++s)
		{
			sOffset = s * sim->nAlleles;
			for (int a = 0; a < sim->nAlleles; ++a)
			{
				int b = sim->alleles2acgt[a];

				if (-1 == b)
				{
					continue;
				}
				sim->fmt_ad_arr[sOffset + a] = acgt_fmt_ad_arr[sOffset + b];
			}
		}
	}

	// set sorted values for FMT_INFO_R tags
	// 		INFO_AD
	//		INFO_ADF
	//		INFO_ADR
	// if (args->addInfoAD||args->addInfoADF||args->addInfoADR)//TODO
	if (args->addInfoAD)
	{
		for (int a = 0; a < sim->nAlleles; ++a)
		{
			int b = sim->alleles2acgt[a];
			if (-1 == b)
			{
				continue;
			}
			sim->info_ad_arr[a] = acgt_info_ad_arr[b];
		}
	}

	sOffset = 0;
	// set sorted values for FMT_NUMBER_G tags
	// 		GL
	//		GP
	//		PL
	for (int s = 0; s < nSamples; ++s)
	{
		double *sample_gls = sim->gl_vals + s * MAX_NGTS;

		sOffset = s * nGenotypes;
		// nGenotypes ==
		// == sim->current_size_bcf_tag_number[FMT_NUMBER_G]
		// == sim->current_size_bcf_tag_number[bcf_tags[GL].n];

		float *sample_gl_arr = sim->gl_arr + sOffset;

		if (0 == sim->fmt_dp_arr[s])
		{
			for (int j = 0; j < nGenotypes; j++)
			{
				bcf_float_set_missing(sample_gl_arr[j]);
			}
		}
		else
		{
			for (int a1 = 0; a1 < sim->nAlleles; ++a1){
				for(int a2=a1; a2<sim->nAlleles; ++a2){

					int b1= sim->alleles2acgt[a1];
					int b2 = sim->alleles2acgt[a2];

					if (b1 == -1 || b2 == -1)
						continue;

					int newoffset = bcf_alleles2gt(a1, a2);
					int oldoffset = bcf_alleles2gt(b1,b2);

					ASSERT(oldoffset > -1);
					ASSERT(newoffset > -1);
					sample_gl_arr[newoffset] = (float)sample_gls[lut_myGtIdx_to_vcfGtIdx[oldoffset]];
				}
			}
			rescale_likelihood_ratio(sample_gl_arr, nGenotypes);
		}
	}


	// [BEGIN CONSTRUCT QS TAG] ------------------------------------------------------- //
	//
	// example: at site 0
	// sample 1:
	//      sampled base,qual
	//             C,4
	//             T,4
	//             A,1
	//             T,6
	//
	// sample 2:
	//      sampled base,qual
	//             G,7
	//
	// sample 1 sum of quals = 15
	// normalized qsum for A,C,G,T
	// sample1: 1/15,4/15,0,10/15
	// sample2: 0,0,7/7,0
	// qsum: 1/15,4/15,7/7,10/15

	if (1 == args->addQS)
	{
		int a = -1;
		float sum = 0.0;

		// -- calculation of the qsum -- //
		// sum the normalized qsum across all samples
		// to account for differences in coverage
		// modified from: bcftools/bam2bcf.c
		for (int s = 0; s < nSamples; ++s)
		{

			sum = 0.0;
			for (int b = 0; b < 4; ++b)
			{
				sum += samples_acgt_qs[b][s];
			}

			if (0 != sum)
			{
				for (int b = 0; b < 4; ++b)
				{
					a = sim->acgt2alleles[b];
					if (-1 == a)
					{
						continue;
					}

					// per-sample normalization
					// save as ordered (base b --> allele index a)
					sim->qs_arr[a] += (float)(samples_acgt_qs[b][s] / sum);
				}
			}
		}
	}

	// [END CONSTRUCT QS TAG] --------------------------------------------------------- //




	if (1 == args->addPL)
	{

		int x;
		for (int i = 0; i < sim->current_size_bcf_tag_number[bcf_tags[PL].n]; ++i)
		{
			if (bcf_float_is_missing(sim->gl_arr[i]))
			{
				sim->pl_arr[i] = bcf_int32_missing;
			}
			else if (bcf_float_is_vector_end(sim->gl_arr[i]))
			{
				NEVER;
			}
			else
			{
				// TODO alias here too
				((x = lroundf(-10 * sim->gl_arr[i])) > MAXPL) ? sim->pl_arr[i] = MAXPL : sim->pl_arr[i] = x;
				// sim->pl_arr[i] = lroundf(-10 * sim->gl_arr[i]); // when max threshold disabled, vcfgl gives same results as bcftools +tag2tag -- --GL-to-PL
			}
		}
	}


	if (1 == args->addGP)
	{
		for (int s = 0; s < nSamples; ++s)
		{
			sOffset = s * nGenotypes;
			// float *gpp = sim->gp_arr + i * nGenotypes;
			// float *glp = sim->gl_arr + i * nGenotypes;
			float sum = 0;
			for (int j = 0; j < nGenotypes; j++)
			{
				if (bcf_float_is_missing(sim->gl_arr[sOffset + j]))
				{
					bcf_float_set_missing(sim->gp_arr[sOffset + j]);
					continue;
				}
				else if (bcf_float_is_vector_end(sim->gl_arr[sOffset + j]))
				{
					NEVER; // we never expect to truncate the vector and finish early
				}
				else
				{
					sim->gp_arr[sOffset + j] = sim->gl_arr[sOffset + j];
					sim->gp_arr[sOffset + j] = pow(10, sim->gp_arr[sOffset + j]);
					sum += sim->gp_arr[sOffset + j];
				}

			}
			if (sum <= 0)
				continue;
			for (int j = 0; j < nGenotypes; j++)
			{
				if (bcf_float_is_missing(sim->gp_arr[sOffset + j]))
				{
					continue;
				}
				else if (bcf_float_is_vector_end(sim->gp_arr[sOffset + j]))
				{
					NEVER;
				}else{
					sim->gp_arr[sOffset + j] /= sum;
				}
			}
		}
	}


	sim->set_tag_I16();


	// [BEGIN REORDER GENOTYPES] ------------------------------------------------------- //
	// reorder genotypes according to the new alleles order

	int alleles_gts[sim->nHaplotypes];
	for (int i = 0; i < sim->nHaplotypes; ++i)
	{
		alleles_gts[i] = -1;
	}
	for (int s = 0; s < nSamples; s++)
	{
		sOffset = s * SIM_PLOIDY;

		// Assume REF=A ALT=C, thus bcf_gt_allele 0 -> A, 1 -> C
		int a1 = sim->acgt2alleles[bcf_gt_allele(sim->gt_arr[sOffset + 0])];
		int a2 = sim->acgt2alleles[bcf_gt_allele(sim->gt_arr[sOffset + 1])];
		if (-1 == a1 || -1 == a2)
		{
			continue;
		}
		alleles_gts[sOffset + 0] = bcf_gt_phased(a1);
		alleles_gts[sOffset + 1] = bcf_gt_phased(a2);
	}
	bcf_update_genotypes(sim->hdr, rec, alleles_gts, sim->nHaplotypes);
	// [END REORDER GENOTYPES] --------------------------------------------------------- //


	sim->add_tags();
	// TODO with trimming enabled, if the allele in the real genotype is not observed, and gt tag is not set to drop, give error and exit

	return (0);
}

int simulate_record_true_values(simRecord *sim)
{

	bcf1_t *rec = sim->rec;

	const int nSamples = sim->nSamples;
	const int nGenotypes = sim->nGenotypes;

	sim->current_size_bcf_tag_number[FMT_NUMBER_G] = sim->nSamples * nGenotypes;

	int32_t ngt_arr = 0;
	int ngt = bcf_get_genotypes(sim->hdr, rec, &sim->gt_arr, &ngt_arr);
	if (ngt <= 0)
	{
		ERROR("Could not find GT tag.");
	}
	if (2 != ngt / nSamples)
	{
		ERROR("Ploidy %d is not supported.\n", ngt / nSamples);
	}

	sim->reset_rec_objects();

	for (int sample_i = 0; sample_i < nSamples; sample_i++)
	{

		int32_t *gt_ptr = sim->gt_arr + sample_i * SIM_PLOIDY;

		// bcf_gt_allele(gt_ptr[i])

		int bin_gts[2] = {0};

		// binary input genotypes from simulated input
		for (int i = 0; i < SIM_PLOIDY; i++)
		{
			bin_gts[i] = bcf_gt_allele(gt_ptr[i]);
			// use bit shifting to check if bin_gt[n] is 0 or 1
			if ((bin_gts[i] >> 1) != 0)
			{
				fprintf(stderr, "ERROR:\n\nbin_gts[%d]: Genotype %d not supported.\n", i, bin_gts[i]);
				exit(1);
			}
		}

		// most likely value should be at the genotype combination for true genotype
		int true_gt_idx = bin_gts[0] + bin_gts[1];
		for (int i = 0; i < nGenotypes; ++i)
		{
			if (i == true_gt_idx)
			{
				sim->gl_arr[sample_i * nGenotypes + i] = MAXGL;

				if (1 == args->addPL)
				{
					sim->pl_arr[sample_i * nGenotypes + i] = MAXPL;
				}
				if (1 == args->addGP)
				{
					sim->gp_arr[sample_i * nGenotypes + i] = MAXGP;
				}
			}
			// else case is already handled by setting all values to MINXX during bcf_tag_reset
		}

	} // end sample loop

	// set REF and ALT alleles
	ASSERT(0 == (bcf_update_alleles_str(sim->hdr, rec, "A,C")));

	if (1 == args->addGL)
	{
		ASSERT(0 == (bcf_update_format_float(sim->hdr, rec, "GL", sim->gl_arr, sim->current_size_bcf_tag_number[bcf_tags[GL].n])));
	}

	if (1 == args->addGP)
	{
		ASSERT(0 == (bcf_update_format_float(sim->hdr, rec, "GP", sim->gp_arr, sim->current_size_bcf_tag_number[bcf_tags[GP].n])));
	}

	if (1 == args->addPL)
	{
		ASSERT(0 == (bcf_update_format_int32(sim->hdr, rec, "PL", sim->pl_arr, sim->current_size_bcf_tag_number[bcf_tags[PL].n])));
	}

	return (0);
}

int main(int argc, char **argv)
{

	args = args_get(--argc, ++argv);

	char *in_fn = args->in_fn;
	const int pos0 = args->pos0;

	FILE *arg_ff = openFILE(args->out_fnp, ".arg");
	ASSERT(NULL != arg_ff);
	
	// fprintf(stderr, "vcfgl [version: %s] [build: %s %s] [htslib: %s]\n", VCFGL_VERSION, __DATE__, __TIME__, hts_version());
	// fprintf(arg_ff, "vcfgl [version: %s] [build: %s %s] [htslib: %s]\n", VCFGL_VERSION, __DATE__, __TIME__, hts_version());


	fprintf(stderr, "\n%s", args->command);
	fprintf(arg_ff, "\n%s", args->command);

	fprintf(stderr, "\nReading file:\t\"%s\"\n", in_fn);
	htsFile *in_ff = hts_open(in_fn, "r");
	ASSERT(NULL != in_ff);

	bcf_hdr_t *in_hdr = bcf_hdr_read(in_ff);
	simRecord *sim = new simRecord(in_hdr);

	fprintf(stderr, "\n\t-> Opening output file for writing: %s\n", args->out_fn);
	htsFile* out_ff = hts_open(args->out_fn, args->output_mode_str);

	// create multithreaded pool
	htsThreadPool tpool = {NULL, 0};
	if (args->n_threads > 1) {
		tpool.pool = hts_tpool_init(args->n_threads);
		ASSERT(NULL != tpool.pool);
		// add input stream to the pool
		hts_set_opt(in_ff, HTS_OPT_THREAD_POOL, &tpool);
		// add output stream to the pool
		hts_set_opt(out_ff, HTS_OPT_THREAD_POOL, &tpool);
	}

	ASSERT(bcf_hdr_write(out_ff, sim->hdr) == 0);

	fprintf(stderr, "\n%s\n", args->datetime);

	bcf1_t *in_rec = bcf_init();
	int nSites = 0;
	int nSitesSkipped = 0;

	args->error_rate_q = floor(-10 * log10(args->error_rate));

	// /BEGIN/ main sites loop ---------------------------------------------


	//TODO HERE
	if (-999 == args->mps_depth)
	{

		// without error, only 3 possible genotypes
		// input:
		//	REF	ALT
		//	0	1
		//
		// gts:
		// 00 01 11
		sim->nGenotypes = 3;
		sim->current_size_bcf_tag_number[FMT_NUMBER_G] = sim->nSamples * sim->nGenotypes;
		simulate_record = &simulate_record_true_values;
	}
	else
	{
		simulate_record = &simulate_record_values;
	}


	if(1 == args->addI16) {
		get_strand=&sample_strand;
	}else{
		get_strand=&set_strand_to_forward;
	}

	int ret; // return value

	if (0 == args->explode)
	{

		while (0 == bcf_read(in_ff, in_hdr, in_rec))
		{

			sim->rec = bcf_copy(sim->rec, in_rec);

			ret = simulate_record(sim);
			sim->rec->pos += pos0;

			if (ret < 0)
			{
				nSitesSkipped++;
				continue;
			}

			ASSERT(0 == bcf_write(out_ff, sim->hdr, sim->rec));
			nSites++;
		}

	}
	else if (1 == args->explode)
	{
		// EXPLODE 1 -------------------------------------------------------- //
		// simulate by enumerating the invariable sites that are not present in the input vcf

		// read the first record
		ASSERT(0 == bcf_read(in_ff, in_hdr, in_rec));

		// prepare blank record template
		bcf1_t *blank_rec = bcf_init();
		blank_rec = bcf_copy(blank_rec, in_rec);

		int32_t *tmp_gt_arr = NULL;
		int32_t tmp_n_gt_arr = 0;
		int tmp_ngt = bcf_get_genotypes(in_hdr, blank_rec, &tmp_gt_arr, &tmp_n_gt_arr);
		if (tmp_ngt <= 0)
		{
			ERROR("Could not find GT tag.");
		}

		for (int h = 0; h < sim->nHaplotypes; ++h)
		{
			tmp_gt_arr[h] = BCF_GT_PHASED_0;
		}

		do
		{

			while (1)
			{ // this block should run for every site with missing/nodata
				if (nSites == in_rec->pos + pos0)
				{
					break;
				}

				blank_rec->pos = nSites;
				ASSERT(0 == (bcf_update_genotypes(in_hdr, blank_rec, tmp_gt_arr, sim->nHaplotypes)));

				sim->rec = bcf_copy(sim->rec, blank_rec);

				ret = simulate_record(sim);
				if (ret < 0)
				{
					nSitesSkipped++;
					//TODO 
					// nSites++;
				// DEVPRINT("nsites %d",nSites);
						// DEVPRINT("pos %d",in_rec->pos + pos0);
					continue;
				}

				ASSERT(0 == bcf_write(out_ff, sim->hdr, sim->rec));
				nSites++;
			}

				sim->rec = bcf_copy(sim->rec, in_rec);
				ret = simulate_record(sim);
				sim->rec->pos += pos0;
				if (ret < 0)
				{
					nSitesSkipped++;
					continue;
				}
				ASSERT(0 == bcf_write(out_ff, sim->hdr, sim->rec));
				nSites++;

		} while (bcf_read(in_ff, in_hdr, in_rec) == 0);

		bcf_idpair_t *ctg = in_hdr->id[BCF_DT_CTG];
		int contigsize = ctg[in_rec->rid].val->info[0];

		while (nSites < contigsize)
		{
			blank_rec->pos = nSites;
			ASSERT(0 == (bcf_update_genotypes(in_hdr, blank_rec, tmp_gt_arr, sim->nHaplotypes)));
			sim->rec = bcf_copy(sim->rec, blank_rec);
			ret = simulate_record(sim);

			if (ret < 0)
			{
				nSitesSkipped++;
				continue;
			}

			ASSERT(0 == bcf_write(out_ff, sim->hdr, sim->rec));
			nSites++;
		}

		bcf_destroy(blank_rec);
		free(tmp_gt_arr);
		tmp_gt_arr = NULL;
	}

	// /END/ main sites loop -----------------------------------------------

	fprintf(stderr, "Total number of sites simulated: %i\n", nSites);
	fprintf(stderr, "Number of sites skipped: %i\n", nSitesSkipped);

	bcf_hdr_destroy(in_hdr);
	bcf_destroy(in_rec);

	ASSERT(0 == hts_close(in_ff));
	ASSERT(0 == hts_close(out_ff));

	if(NULL!=tpool.pool){
		hts_tpool_destroy(tpool.pool);
	}
	

	ASSERT(0 == fclose(arg_ff));


	delete sim;
	
	args_destroy(args);

	return 0;
}
