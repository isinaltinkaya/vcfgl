/*
 * vcfgl: gl simulator for vcf/bcf files
 *
 * isinaltinkaya
 *
 *
 */

#include <stdio.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <inttypes.h>
#include <math.h>
#include <time.h>

#include "shared.h"
#include "random_generator.h"
#include "estimator.h"
#include "io.h"
#include "bcf_utils.h"

argStruct *args;

/// @brief sample_strand - sample a strand
/// @return		0: forward, 1: reverse
inline int sample_strand()
{
	return (drand48() < 0.5 ? 0 : 1);
}

/// @brief sample_tail_distance - sample a tail length for I16 tag
/// @return		int length
inline int sample_tail_distance()
{
	int i;
	if ((i = sample_uniform_from_range(1, 50)) > CAP_DIST)
	{
		return (CAP_DIST);
	}
	return (i);
}

int create_blank_record(sim_rec *sim, bcf1_t *blk, bcf1_t *unmod, bcf_hdr_t *hdr)
{

	int32_t ngt_arr = 0;
	int ngt = bcf_get_genotypes(hdr, blk, &sim->gt_arr, &ngt_arr);
	if (ngt <= 0)
	{
		fprintf(stderr, "\nGT not present\n");
	}

	int32_t *tmpia = sim->gt_arr;
	for (int i = 0; i < bcf_hdr_nsamples(hdr); i++)
	{
		tmpia[2 * i + 0] = bcf_gt_phased(0);
		tmpia[2 * i + 1] = bcf_gt_phased(0);
	}
	ASSERT(0 == (bcf_update_genotypes(hdr, blk, tmpia, bcf_hdr_nsamples(hdr) * 2)));

	return 0;
}

int simulate_record_values(sim_rec *sim, FILE *out_baseCounts_ff)
{

	// -------------------------------------------------------------------------
	// reset tag arrays with fixed size
	bcf_tag_reset<int32_t>(sim->dp_arr, DP, 0);
	if (1 == args->addQS)
	{
		bcf_tag_reset<float>(sim->qs_arr, QS, 0);
	}

	// -------------------------------------------------------------------------
	bcf_tag_reset<float>(sim->gl_arr, GL, -0.0);

	if (1 == args->addGP)
	{
		bcf_tag_reset<float>(sim->gp_arr, GP, -0.0);
	}

	if (1 == args->addPL)
	{
		bcf_tag_reset<int32_t>(sim->pl_arr, PL, 0);
	}
	// -------------------------------------------------------------------------

	// -------------------------------------------------------------------------

	bcf1_t *rec = NULL;
	(sim->blank_rec == NULL) ? rec = sim->rec : rec = sim->blank_rec;

	kstring_t *alleles_str = kbuf_init();

	const int nSamples = sim->nSamples;

	double *gls = sim->gl_vals;

	int nGenotypes = -1;

	// TODO delme
	for (int i = 0; i < nSamples * MAX_NGTS; ++i)
	{
		sim->gl_vals[i] = -0.0;
	}

	int n_sim_reads = -1;

	int qScore = -1;

	int which_strand, which_haplo = -1;
	int refnref = -1;
	int e_base, base = -1;
	double e = -1.0; // simulated error rate instance

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

	// refnref_n_q13_bases	number of bases where quality score >= 13 summed across all individuals
	//
	// \def refnref_n_q13_bases[2][2]
	// 		refnref_n_q13_bases[reference|non-reference][forward-strand|reverse-strand]
	//
	// e.g.
	// 		refnref_n_q13_bases[0][0] = #reference bases on forward-strand with q >= 13
	int refnref_n_q13_bases[2][2] = {{0, 0}, {0, 0}};

	// refnref_sum_taildist	sum of tail lengths
	// e.g.
	// 		refnref_sum_taildist[1]	= sum of tail lengths for non-ref bases
	//
	int refnref_sum_taildist[2] = {0, 0};
	int refnref_sum_taildist_sq[2] = {0, 0};

	int refnref_sum_qs[2] = {0, 0};
	int refnref_sum_qs_sq[2] = {0, 0};

	// \def acgt_counts[4]	base counts
	int acgt_counts[4] = {0};

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

	// acgt2alleles[i] gives the index of the allele in ordered bcf REF+ALT alleles
	// i		index of the allele in {A,C,G,T}
	// return	index of the allele in ordered bcf REF+ALT alleles
	// e.g.
	// allele			A,C,G,T
	// acgt_index i 	0,1,2,3
	// ordered_index 	3,1,0,2 for a site where REF=T ALT=C,A,G
	// for qs reordering
	int acgt2alleles[5]; // alleles: ref, alt, alt2, alt3...
	acgt2alleles[0] = acgt2alleles[1] = acgt2alleles[2] = acgt2alleles[3] = acgt2alleles[4] = -1;

	// A,C,G,T,N
	int seens[5] = {0, 1, 2, 3, 4};

	int tail_dist = 0;

	for (int sample_i = 0; sample_i < nSamples; sample_i++)
	{

		if (args->printBaseCounts == 1)
		{
			acgt_counts[0] = acgt_counts[1] = acgt_counts[2] = acgt_counts[3] = 0;
		}

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
			sim->dp_arr[sample_i] = 0;
		}
		else
		{

			int32_t *gt_ptr = sim->gt_arr + sample_i * SIM_PLOIDY;

			double *sample_gls = gls + sample_i * MAX_NGTS;

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

			int sim_strands[n_sim_reads] = {-1};
			// int sim_quals[n_sim_reads]={-1};

			for (int i = 0; i < n_sim_reads; i++)
			{

				(drand48() < 0.5) ? which_haplo = 0 : which_haplo = 1;

				e = drand48();

				e_base = base = bin_gts[which_haplo];
				if (e < args->errate)
				{
					do
					{
						e_base = floor(4 * drand48());
					} while (e_base == base);
					base = e_base;
				}

				// in expected input file, ref is always 0 and non-ref is always 1
				// refnref == 0 -> we have ref
				// refnref == 1 -> we have nonref
				(0 == base) ? refnref = 0 : refnref = 1;

				gl_log10(base, args->errate, sample_gls);

				acgt_counts[base]++;

				qScore = floor(-10 * log10(e));
				if (qScore > CAP_BASEQ)
					qScore = CAP_BASEQ;

				(1 == args->addI16) ? which_strand = sample_strand() : which_strand = 1;
				sim_strands[i] = which_strand;

				if (qScore >= 13)
				{
					refnref_n_q13_bases[refnref][which_strand]++;
				}

				refnref_sum_qs[refnref] = qScore;
				refnref_sum_qs_sq[refnref] = qs_to_qs2(qScore);

				tail_dist = sample_tail_distance();
				refnref_sum_taildist[refnref] += tail_dist;
				refnref_sum_taildist_sq[refnref] += tail_dist * tail_dist;

				samples_acgt_qs[base][sample_i] += qScore;
			}

			sim->dp_arr[sample_i] = n_sim_reads;
		}

		if (NULL != out_baseCounts_ff)
		{
			fprintf(out_baseCounts_ff, "%d\t%d\t%d\t%d\t%d\t%d\n", sim->site_i, sample_i, acgt_counts[0], acgt_counts[1], acgt_counts[2], acgt_counts[3]);
		}

	} // end sample loop

	int n_alleles = -1;
	int use_as_alt = -1;

	// [BEGIN QS TAG] ------------------------------------------------------- //
	if (1 == args->addQS)
	{

		// format after seq_nt16_str:
		// 0,1,2,3 for A,C,G,T; 4 otherwise. <0 indel
		// int ori_ref=seq_nt16_int[seq_nt16_table[SIM_acgt2alleles[0]]];
		int ori_ref = 0; // equivalent to above

		int ref4 = 0; // TODO

		int unseen = -1;
		acgt2alleles[0] = ref4;

		float qsum[5] = {0.0};

		// -- calculation of the qsum -- //
		// sum the normalized qsum across all samples
		// to account for differences in coverage
		// modified from: bcftools/bam2bcf.c
		for (int i = 0; i < nSamples; ++i)
		{

			float sum = 0.0;
			for (int j = 0; j < 4; ++j)
			{
				sum += samples_acgt_qs[j][i];
			}

			if (0 != sum)
			{
				for (int j = 0; j < 4; ++j)
				{
					// per-sample normalization
					qsum[j] += (float)samples_acgt_qs[j][i] / sum;
				}
			}
			else
			{
				// WARNING("Sum of the quality scores at site with index %d is %d.",site_i);
			}
		}

		// sort qsum in ascending order (insertion sort)
		float *ptr[5], *tmp;

		for (int i = 0; i < 5; i++)
		{
			ptr[i] = &qsum[i];
		}

		for (int i = 1; i < 4; i++)
		{
			for (int j = i; j > 0 && *ptr[j] < *ptr[j - 1]; j--)
			{
				tmp = ptr[j], ptr[j] = ptr[j - 1], ptr[j - 1] = tmp;
			}
		}

		int i, j = 0;
		// Set the reference allele and alternative allele(s)
		for (i = 3, j = 1; i >= 0; i--) // i: alleles sorted by QS; j, acgt2alleles[j]: output allele ordering
		{
			int ipos = ptr[i] - qsum; // position in sorted qsum array
			ASSERT(ipos > -1 && ipos < 5);
			if (ipos == ref4)
			{
				sim->qs_arr[0] = qsum[ipos]; // REF's qsum
			}
			else
			{
				// NB this will cause the i to be >0 in the below checks
				if (!qsum[ipos])
				{
					// NEVER;
					break; // qsum is 0, this and consequent alleles are not seen in the pileup
				}
				sim->qs_arr[j] = qsum[ipos];
				acgt2alleles[j++] = ipos;
			}
		}
		int ref_base = 0; // TODO delme unnec var & ifs nextline
		if (ref_base >= 0)
		{

			// for SNPs, find the "unseen" base
			if (((ref4 < 4 && j < 4) || (ref4 == 4 && j < 5)) && i >= 0)
			{
				unseen = j, acgt2alleles[j++] = ptr[i] - qsum;
			}
			// observed: T -> n_alleles=2 (T,<*>)
			// observed: A,G -> n_alleles=3 (A,G,<*>)
			// observed: A,C,G -> n_alleles=4 (A,C,G,<*>)
			// observed: A,C,G,T -> n_alleles=4 (A,C,G,T)
			n_alleles = j;
		}
		else
		{
			NEVER;
			// n_alleles = j;
			// // if (n_alleles == 1) return -1; // no reliable supporting read. stop doing anything
			// if (n_alleles == 1) NEVER;
		}
		int has_alt = (n_alleles == 2 && unseen != -1) ? 0 : 1;

		int nseen = 0;	 // #seen bases
		int nunseen = 0; // #unseen bases
		for (int ii = 0; ii < 5; ++ii)
		{
			if ((acgt2alleles[ii] >= 0) && (unseen != ii))
			{
				seens[acgt2alleles[ii]] = -999;
				++nseen;
			}
			else
			{
				nunseen++;
			}
		}

		if (1 == args->trimAlts)
		{
			n_alleles = nseen;
		}
		else if (0 == args->trimAlts)
		{
			n_alleles = 4;
		}

		// @@
		for (int i = 1; i < 5; ++i)
		{
			if (acgt2alleles[i] < 0)
				break;

			if (1 == args->trimAlts)
			{

				if (1 == n_alleles)
				{
					// invar site
					// therefore set ALT to one of REF or ALT
					// e.g.
					// normally in -explode 1 we assume REF=A ALT=C
					// at a site, if only C is observed (therefore 1==n_alleles)
					// 		set ALT to A
					// if only A is observed
					// 		set ALT to C
					// if only G is observed (very unlikely case, only if errate==1
					// 		or very high and by chance we observe only an error)
					// 		set ALT to A (basically first base in REFALT {A,C})

					use_as_alt = 0; // default to first base in REFALT == A

					// exclude the beginning of unseens (so called unseen) from our lookup acgt
					acgt2alleles[unseen] = -1;

					for (int ii = 0; ii < 4; ++ii)
					{ // 4->exclude N
						if (-999 == seens[ii])
						{
							if (ii == 0)
							{
								use_as_alt = 1;
								// DEVPRINT("%d %d %d %d",acgt2alleles[0],acgt2alleles[1],acgt2alleles[2],acgt2alleles[3]);
								// int old=acgt2alleles[use_as_alt];
								acgt2alleles[use_as_alt] = use_as_alt;
								// DEVPRINTX("ii=%d use_as_alt=%d acgt2alleles[use_as_alt] old =%d new=%d",ii,use_as_alt,old,acgt2alleles[use_as_alt]);
							}
							else if (ii == 1)
							{
								use_as_alt = 0;
								acgt2alleles[use_as_alt] = use_as_alt;
							}
							else
							{
								NEVER;
							}
						}
					}

				}
			}
		}
		//
		//
		// @@







		int nals = 1;
		kputc("ACGTN"[ori_ref], alleles_str);
		for (int i = 1; i < 5; ++i)
		{
			if (acgt2alleles[i] < 0)
				break;

			if (1 == args->trimAlts)
			{

				if (1 == n_alleles)
				{
					n_alleles = 2; // n_alleles is now 2 since we set ALT above
					kputc(',', alleles_str);
					kputc("ACGT"[use_as_alt], alleles_str);
				}
				else
				{

					// trim enabled but all bases were observed
					if(unseen==-1){
							kputc(',', alleles_str);
							kputc("ACGT"[acgt2alleles[i]], alleles_str);
					}else{

						if (i != 4 && acgt2alleles[i + 1] >= 0)
						{
							kputc(',', alleles_str);
						}


						if (unseen == i)
						{
							// kputs("<*>", alleles_str);
						}
						else
						{
							ASSERT(unseen!=i);
							kputc("ACGT"[acgt2alleles[i]], alleles_str);
						}
					}
					//@@
					// kputc("ACGT"[acgt2alleles[i]], alleles_str);
					//@@
				}
			}
			else if (0 == args->trimAlts)
			{

				kputc(',', alleles_str);

				if (unseen == i)
				{

					int x = nunseen;
					for (int ii = 0; ii < 4; ++ii)
					{ // 4->exclude N
						if (-999 == seens[ii])
						{
							// DEVPRINT("seen %d -> %c",ii,"ACGTN"[ii]);
						}
						else
						{
							// DEVPRINT("unseen %d -> %c",ii,"ACGTN"[ii]);

							if (x > 0 && x != nunseen)
							{
								kputc(',', alleles_str);
							}
							x--;

							kputc("ACGT"[ii], alleles_str);
						}
					}
				}
				else
				{
					kputc("ACGT"[acgt2alleles[i]], alleles_str);
				}
			}
			else
			{
				NEVER;
			}
		}
	}

	// [END QS TAG] --------------------------------------------------------- //

	if (1 == args->trimAlts)
	{
		nGenotypes = nAlleles_to_nGenotypes(n_alleles);
	}
	else
	{
		nGenotypes = MAX_NGTS;
	}

	const int sizeToUse = nSamples * nGenotypes;

	// set REF and ALT alleles
	if (0 == alleles_str->l)
	{
		ASSERT(1 != args->addQS);
		kputs("A,C,G,T", alleles_str);
	}
	// ASSERT(alleles_str->l>1);
	ASSERT(0 == (bcf_update_alleles_str(sim->hdr, rec, alleles_str->s)));

	ASSERT(0 == (bcf_update_format_int32(sim->hdr, rec, "DP", sim->dp_arr, bcf_tags[DP].n)));

	// [BEGIN GL TAG] ------------------------------------------------------- //
	for (int sample_i = 0; sample_i < nSamples; ++sample_i)
	{

		double *sample_gls = gls + sample_i * MAX_NGTS;
		// rescale_likelihood_ratio(sample_gls);

		float *sample_gl_arr = sim->gl_arr + sample_i * nGenotypes;

		if (0 == sim->dp_arr[sample_i])
		{
			if (0 == args->addQS)
			{
				ASSERT(nGenotypes == MAX_NGTS);
			}
			for (int j = 0; j < nGenotypes; j++)
			{
				bcf_float_set_missing(sim->gl_arr[sample_i * nGenotypes + j]);
			}
		}
		else
		{
			if (0 == args->trimAlts)
			{
				ASSERT(nGenotypes == MAX_NGTS);
				for (int j = 0; j < nGenotypes; j++)
				{
					sample_gl_arr[j] = (float)sample_gls[lut_myGtIdx_to_vcfGtIdx[j]];
				}
			}
			else
			{

				for (int a1 = 0; a1 < 4; ++a1)
				{
					for (int a2 = a1; a2 < 4; ++a2)
					{

						int aa1 = acgt2alleles[a1];
						int aa2 = acgt2alleles[a2];

						if (aa1 == -1 || aa2 == -1)
							continue;
						// TODO

						int newoffset = bcf_alleles2gt(a1, a2);
						int oldoffset = bcf_alleles2gt(aa1, aa2);

						ASSERT(oldoffset > -1);
						ASSERT(newoffset > -1);
						sample_gl_arr[newoffset] = (float)sample_gls[lut_myGtIdx_to_vcfGtIdx[oldoffset]];
					}
				}
			}
			rescale_likelihood_ratio(sample_gl_arr,nGenotypes);
		}
	}
	gls = NULL;

	kbuf_destroy(alleles_str);

	ASSERT(0 == (bcf_update_format_float(sim->hdr, rec, "GL", sim->gl_arr, sizeToUse)));
	// [END GL TAG] --------------------------------------------------------- //

	if (1 == args->addPL)
	{

		int x;
		for (int i = 0; i < bcf_tags[PL].n; ++i)
		{
			if (bcf_float_is_missing(sim->gl_arr[i]))
			{
				sim->pl_arr[i] = bcf_int32_missing;
			}
			else if (bcf_float_is_vector_end(sim->gl_arr[i]))
			{
				NEVER;
				// sim->pl_arr[i] = bcf_int32_vector_end;
			}
			else
			{
				// ((x = lroundf(-10 * sim->gl_arr[i])) > MAXPL) ? sim->pl_arr[i] = MAXPL : sim->pl_arr[i] = x;
				x = lroundf(-10 * sim->gl_arr[i]);
				// (x > MAXPL) ? sim->pl_arr[i]=MAXPL:sim->pl_arr[i]=x;
				if(x>MAXPL){

					sim->pl_arr[i]=MAXPL;
				}else{
					sim->pl_arr[i]=x;
				}
			}
		}
		ASSERT(0 == (bcf_update_format_int32(sim->hdr, rec, "PL", sim->pl_arr, sizeToUse)));
	}

	if (1 == args->addGP)
	{
		for (int i = 0; i < nSamples; i++)
		{
			float *gpp = sim->gp_arr + i * nGenotypes;
			float *glp = sim->gl_arr + i * nGenotypes;
			float sum = 0;
			for (int j = 0; j < nGenotypes; j++)
			{
				if (bcf_float_is_vector_end(glp[j]))
				{
					NEVER; // we never expect to truncate the vector and finish early
				}
				if (bcf_float_is_missing(glp[j]))
				{
					bcf_float_set_missing(gpp[j]);
					continue;
				}
				else
				{
					gpp[j] = glp[j];
				}
				gpp[j] = pow(10, gpp[j]);
				sum += gpp[j];
			}
			if (sum <= 0)
				continue;
			for (int j = 0; j < nGenotypes; j++)
			{
				if (bcf_float_is_missing(gpp[j]))
				{
					continue;
				}
				if (bcf_float_is_vector_end(gpp[j]))
				{
					NEVER;
					// break;
				}

				gpp[j] /= sum;
			}
		}
		ASSERT(0 == (bcf_update_format_float(sim->hdr, rec, "GP", sim->gp_arr, sizeToUse)));
	}

	if (1 == args->addI16)
	{

		// refnref_n_q13_bases[][]
		// [ref|nonref][forward|reverse]

		// 1   #reference Q13 bases on the forward strand
		sim->i16_arr[0] = refnref_n_q13_bases[0][0];

		// 2   #reference Q13 bases on the reverse strand
		sim->i16_arr[1] = refnref_n_q13_bases[0][1];

		// 3   #non-ref Q13 bases on the forward strand
		sim->i16_arr[2] = refnref_n_q13_bases[1][0];

		// 4   #non-ref Q13 bases on the reverse strand
		sim->i16_arr[3] = refnref_n_q13_bases[1][1];

		// 5   sum of reference base qualities
		sim->i16_arr[4] = refnref_sum_qs[0];

		// 6   sum of squares of reference base qualities
		sim->i16_arr[5] = refnref_sum_qs_sq[0];

		// 7   sum of non-ref base qualities
		sim->i16_arr[6] = refnref_sum_qs[1];

		// 8   sum of squares of non-ref base qualities
		sim->i16_arr[7] = refnref_sum_qs_sq[1];

		// 9   sum of ref mapping qualities
		sim->i16_arr[8] = 0;

		// 10  sum of squares of ref mapping qualities
		sim->i16_arr[9] = 0;

		// 11  sum of non-ref mapping qualities
		sim->i16_arr[10] = 0;

		// 12  sum of squares of non-ref mapping qualities
		sim->i16_arr[11] = 0;

		// 13  sum of tail distance for ref bases
		sim->i16_arr[12] = refnref_sum_taildist[0];

		// 14  sum of squares of tail distance for ref bases
		sim->i16_arr[13] = refnref_sum_taildist_sq[0];

		// 15  sum of tail distance for non-ref bases
		sim->i16_arr[14] = refnref_sum_taildist[1];

		// 16  sum of squares of tail distance for non-ref
		sim->i16_arr[15] = refnref_sum_taildist_sq[1];

		ASSERT(0 == (bcf_update_info_float(sim->hdr, rec, "I16", sim->i16_arr, 16)));
	}

	// [BEGIN ADD QS TAG] ------------------------------------------------------- //
	if (1 == args->addQS)
	{
		ASSERT(0 == (bcf_update_info_float(sim->hdr, rec, "QS", sim->qs_arr, bcf_tags[QS].n)));
	}
	// [END ADD QS TAG] --------------------------------------------------------- //

	return (0);
}

int simulate_record_true_values(sim_rec *sim)
{

	bcf1_t *rec = NULL;
	(sim->blank_rec == NULL) ? rec = sim->rec : rec = sim->blank_rec;

	const int nSamples = sim->nSamples;

	// without error, only 3 possible genotypes
	// input:
	//	REF	ALT
	//	0	1
	//
	// gts:
	// 00 01 11
	const int nGenotypes = 3;
	const int sizeToUse = nSamples * nGenotypes;

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

	bcf_tag_reset<float>(sim->gl_arr, GL, -0.0);

	if (1 == args->addGP)
	{
		bcf_tag_reset<float>(sim->gp_arr, GP, -0.0);
	}

	if (1 == args->addPL)
	{
		bcf_tag_reset<int32_t>(sim->pl_arr, PL, 0);
	}

	for (int sample_i = 0; sample_i < nSamples; sample_i++)
	{

		if (args->printBaseCounts == 1)
		{
			NEVER;
		}

		int32_t *gt_ptr = sim->gt_arr + sample_i * SIM_PLOIDY;

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
			else
			{
				sim->gl_arr[sample_i * nGenotypes + i] = MINGL;

				if (1 == args->addPL)
				{
					sim->pl_arr[sample_i * nGenotypes + i] = MINPL;
				}
				if (1 == args->addGP)
				{
					sim->gp_arr[sample_i * nGenotypes + i] = MINGP;
				}
			}
		}

	} // end sample loop

	// set REF and ALT alleles
	ASSERT(0 == (bcf_update_alleles_str(sim->hdr, rec, "A,C")));

	ASSERT(0 == (bcf_update_format_float(sim->hdr, rec, "GL", sim->gl_arr, sizeToUse)));

	if (1 == args->addGP)
	{
		ASSERT(0 == (bcf_update_format_float(sim->hdr, rec, "GP", sim->gp_arr, sizeToUse)));
	}

	if (1 == args->addPL)
	{
		ASSERT(0 == (bcf_update_format_int32(sim->hdr, rec, "PL", sim->pl_arr, sizeToUse)));
	}

	if (1 == args->addI16)
	{
		NEVER;
	}

	if (1 == args->addQS)
	{
		NEVER;
	}

	return (0);
}

int simulate_record(sim_rec *sim, FILE *out_baseCounts_ff)
{

	if (-999 == args->mps_depth)
	{
		ASSERT(NULL == sim->mps_depths);
		ASSERT(NULL == out_baseCounts_ff);
		return simulate_record_true_values(sim);
	}

	return simulate_record_values(sim, out_baseCounts_ff);
}

int main(int argc, char **argv)
{

	args = args_get(--argc, ++argv);

	if (args != NULL)
	{

		char *in_fn = args->in_fn;
		char *out_fp = args->out_fp;
		char *in_mps_depths = args->in_mps_depths;
		const int pos0 = args->pos0;

		FILE *arg_ff = openFILE(out_fp, ".arg");

		fprintf(stderr, "\n%s", args->command);
		fprintf(arg_ff, "\n%s", args->command);

		FILE *out_baseCounts_ff = NULL;
		if (args->printBaseCounts == 1)
		{
			out_baseCounts_ff = openFILE(out_fp, ".baseCounts.tsv");
			fprintf(out_baseCounts_ff, "site\tind\tA\tC\tG\tT\n");
		}

		vcfFile *in_ff = bcf_open(in_fn, "r");
		vcfFile *out_ff = NULL;

		char *OUT_EXT = NULL;
		char *out_fn = NULL;
		switch (*args->output_mode)
		{
		case 'v':
			fprintf(stderr, "\nOutput is VCF file\n");
			OUT_EXT = strdup(".vcf");
			out_fn = (char *)malloc(strlen(out_fp) + strlen(OUT_EXT) + 1);
			strcpy(out_fn, out_fp);
			strcat(out_fn, OUT_EXT);
			fprintf(stderr, "\n\t-> Opening output file for writing: %s\n", out_fn);
			out_ff = bcf_open(out_fn, "w");
			break;
		case 'b':
			fprintf(stderr, "\nOutput is BCF file\n");
			OUT_EXT = strdup(".bcf");
			out_fn = (char *)malloc(strlen(out_fp) + strlen(OUT_EXT) + 1);
			strcpy(out_fn, out_fp);
			strcat(out_fn, OUT_EXT);
			fprintf(stderr, "\n\t-> Opening output file for writing: %s\n", out_fn);
			out_ff = bcf_open(out_fn, "wb");
			break;
		case 'z':
			fprintf(stderr, "\nOutput is compressed VCF file\n");
			OUT_EXT = strdup(".vcf.gz");
			out_fn = (char *)malloc(strlen(out_fp) + strlen(OUT_EXT) + 1);
			strcpy(out_fn, out_fp);
			strcat(out_fn, OUT_EXT);
			fprintf(stderr, "\n\t-> Opening output file for writing: %s\n", out_fn);
			out_ff = bcf_open(out_fn, "wz");
			break;
		case 'u':
			fprintf(stderr, "\nOutput is uncompressed BCF file\n");
			OUT_EXT = strdup(".bcf");
			out_fn = (char *)malloc(strlen(out_fp) + strlen(OUT_EXT) + 1);
			strcpy(out_fn, out_fp);
			strcat(out_fn, OUT_EXT);
			fprintf(stderr, "\n\t-> Opening output file for writing: %s\n", out_fn);
			out_ff = bcf_open(out_fn, "wbu");
			break;
		}
		free(OUT_EXT);

		if (in_ff == NULL)
		{
			return 1;
		}

		if (bcf == 0)
		{
			return 1;
		}

		bcf_hdr_t *in_hdr = bcf_hdr_read(in_ff);
		const int nSamples = bcf_hdr_nsamples(in_hdr);

		sim_rec *sim = new sim_rec(in_hdr);
		ASSERT(bcf_hdr_write(out_ff, sim->hdr) == 0);

		fprintf(stderr, "\n%s\n", args->datetime);

		bcf1_t *bcf = bcf_init();
		int nSites = 0;

		fprintf(stderr, "\nReading file:\t\"%s\"\n", in_fn);
		fprintf(stderr, "Number of samples: %i\n", nSamples);
		fprintf(stderr, "Number of contigs: %d\n", in_hdr->n[BCF_DT_CTG]);

		bcf1_t *blank = bcf_init();
		bcf1_t *rec = bcf_init();

		// /BEGIN/ main sites loop ---------------------------------------------

		// seperate the first read loop from others to avoid having pos==-1 check in each loop
		ASSERT(0 == bcf_read(in_ff, in_hdr, bcf));
		if (bcf->pos == -1)
		{
			if (pos0 == 0)
			{
				fprintf(stderr, "\n[ERROR]: Input file coordinates start from 0; but -pos0 is not set to 1. Please run again with -pos0 1.\n\n");
				exit(1);
			}
		}

		if (0 == args->explode)
		{

			do
			{
				rec = bcf_copy(rec, bcf);
				sim->rec = rec;

				sim->site_i = nSites;
				ASSERT(0 == simulate_record(sim, out_baseCounts_ff));
				sim->rec->pos += pos0;
				ASSERT(0 == bcf_write(out_ff, sim->hdr, sim->rec));
				nSites++;

			} while (bcf_read(in_ff, in_hdr, bcf) == 0);

			// EXPLODE -------------------------------------------------------------
			// simulate by enumerating the invariable sites that are not present in the input vcf
		}
		else if (1 == args->explode)
		{

			do
			{

				rec = bcf_copy(rec, bcf);
				sim->rec = rec;

				blank = bcf_copy(blank, sim->rec);

				while (1)
				{ // this block should run for every site with missing/nodata
				  //   fprintf(stderr,"\t\t-> out_bcf_recpos: %d nSites: %d\n",sim->rec->pos,nSites);
					if (nSites == sim->rec->pos + pos0)
					{
						// fprintf(stderr,"now breaking\n");
						break;
					}

					create_blank_record(sim, blank, sim->rec, in_hdr);
					sim->blank_rec = blank;
					sim->site_i = nSites;
					ASSERT(0 == simulate_record(sim, out_baseCounts_ff));
					blank->pos = nSites;
					// fprintf(stderr,"blank->pos: %d\n",blank->pos+pos0);
					if (bcf_write(out_ff, sim->hdr, blank) != 0)
					{
						fprintf(stderr, "Error: Failed to write\n");
						exit(1);
					}
					sim->blank_rec = NULL;
					nSites++;
				}
				// fprintf(stderr,"After loop that fills in missing data will print out: %d\n",sim->rec->pos+pos0+1);
				sim->site_i = nSites;
				ASSERT(0 == simulate_record(sim, out_baseCounts_ff));
				sim->rec->pos += pos0;
				if (bcf_write(out_ff, sim->hdr, sim->rec) != 0)
				{
					fprintf(stderr, "Error: Failed to write\n");
					exit(1);
				}
				nSites++;

			} while (bcf_read(in_ff, in_hdr, bcf) == 0);

			bcf_idpair_t *ctg = in_hdr->id[BCF_DT_CTG];
			int contigsize = ctg[sim->rec->rid].val->info[0];

			while (nSites < contigsize)
			{
				create_blank_record(sim, blank, sim->rec, in_hdr);
				sim->blank_rec = blank;
				sim->site_i = nSites;
				ASSERT(0 == simulate_record(sim, out_baseCounts_ff));
				blank->pos = nSites;
				if (bcf_write(out_ff, sim->hdr, blank) != 0)
				{
					fprintf(stderr, "Error: Failed to write\n");
					exit(1);
				}
				sim->blank_rec = NULL;
				nSites++;
			}
		}
		else
		{
			NEVER;
		}

		// /END/ main sites loop -----------------------------------------------

		fprintf(stderr, "Total number of sites: %i\n", nSites);

		bcf_hdr_destroy(in_hdr);
		bcf_destroy(bcf);
		bcf_destroy(blank);
		bcf_destroy(rec);

		ASSERT(0 == bcf_close(in_ff));
		ASSERT(0 == bcf_close(out_ff));
		ASSERT(0 == fclose(arg_ff));

		if (args->printBaseCounts == 1)
		{
			ASSERT(0 == fclose(out_baseCounts_ff));
			fprintf(stderr, "\nDumping baseCounts file to %s.baseCounts.tsv\n", out_fp);
		}

		free(out_fn);
		out_fn = NULL;

		delete sim;

		args_destroy(args);
	}

	return 0;
}
