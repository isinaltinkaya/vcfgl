/*
 * vcfgl: gl simulator for vcf/bcf files
 *
 * isinaltinkaya
 *
 */

#include <htslib/thread_pool.h> // htsThreadPool
#include <time.h> // clock

#include "shared.h"
#include "io.h"
#include "random_generator.h"
#include "bcf_utils.h"
#include "gl_methods.h"

argStruct* args;

const char* nonref_str;


int rec_alleles[5] = { -1, -1, -1, -1, -1 };
int* true_gts_acgt_int = NULL;
int* true_gts_alleles_idx = NULL;
int* n_sim_reads_arr = NULL;

// return value < 0		 skip the site
// returns negative value only if program will skip the position given arg values
// skip position reasons:
// -1   all input true genotypes at position are homozygous ref (iff PROGRAM_WILL_SKIP_INPUT_HOMOREFGT_SITES)
// -2   all input true genotypes at position are homozygous alt (iff PROGRAM_WILL_SKIP_INPUT_HOMOALTGT_SITES)
/// @note skips the site in truth VCF as well
inline int check_rec_alleles(simRecord* sim) {

    // -------------------------------------------- //
    // read genotypes
    int32_t ngt_arr = 0;
    int ngt = 0;
    ngt = bcf_get_genotypes(sim->hdr, sim->rec, &sim->gt_arr, &ngt_arr);
    if (ngt <= 0) {
        ERROR("Could not find GT tag at position %ld.", sim->rec->pos + 1);
    }

    const int n_alleles = sim->rec->n_allele;
    if (n_alleles > 5) {
        ERROR("Multiallelic sites with more than 4 alleles are not supported. If this is a feature you need, please contact the developers or create a feature request on the GitHub page.");
    }

    const int nSamples = sim->nSamples;



    // -> clear
    // rec_alleles[0] = -1;
    // rec_alleles[1] = -1;
    // rec_alleles[2] = -1;
    // rec_alleles[3] = -1;
    // rec_alleles[4] = -1;


    char* allele = NULL;
    int x;
    for (int i = 0;i < n_alleles;++i) {

        allele = sim->rec->d.allele[i];

        if (allele[0] == '<') {
            ERROR("Allele '%s' at position %ld is not a valid base.", allele, sim->rec->pos + 1);
        }
        if (allele[1] != '\0') {
            ERROR("Allele '%s' at position %ld is not a single base. If this is on purpose and this is a feature you want to use, please contact the developers or create a feature request on the GitHub page.", allele, sim->rec->pos + 1);
        }


        if (PROGRAM_WILL_USE_ACGT_GTSOURCE) {
            if (-1 == (rec_alleles[i] = ACGT_CHAR2INT(allele[0]))) {
                ERROR("Allele '%s' at position %ld is not a valid base.", allele, sim->rec->pos + 1);
            }

        } else if (PROGRAM_WILL_USE_BINARY_GTSOURCE) {
            // assume: input vcf always have REF=0 ALT=1 
            // #CHROM	POS	ID	REF	ALT	QUAL    FILTER  INFO    FORMAT  ind
            // chr	    1   .   0   1   .       PASS    .       GT      0|1
            // then set REF to A and ALT to C:
            // chr	    1   .   A   C   .       PASS    .       GT      0|1

            // make sure REF and ALT are binary
            x = allele[0] - '0';
            if (x != 0 && x != 1) {
                ERROR("[--source %d] Found allele '%s' at position %ld. Only 0 and 1 are allowed when using binary GT source. Please change the source or check your input file.", args->gtSource, allele, sim->rec->pos + 1);
            }
            rec_alleles[i] = x;
        }

    }

    if (PROGRAM_WILL_USE_BINARY_GTSOURCE) {

        if (n_alleles > 2) {
            ERROR("Multiallelic sites are not supported when using binary GT source. If this is a feature you need, please contact the developers or create a feature request on the GitHub page.");
        }

        ASSERT(0 == (bcf_update_alleles_str(sim->hdr, sim->rec, "A,C")));
    }

    int allelesum = 0;
    int a;
    for (int i = 0;i < nSamples * SIM_PLOIDY;++i) {
        true_gts_acgt_int[i] = -1; // clear

        a = bcf_gt_allele(sim->gt_arr[i]);
        true_gts_alleles_idx[i] = a;
        allelesum += a;

        true_gts_acgt_int[i] = ACGT_CHAR2INT(sim->rec->d.allele[a][0]);
    }

    if (PROGRAM_WILL_SKIP_INPUT_HOMOREFGT_SITES && (0 == allelesum)) {
        return(-1);
    }

    // for each alt, check if all inds are homozygous alt
    if (PROGRAM_WILL_SKIP_INPUT_HOMOALTGT_SITES) {
        for (a = 1; a < n_alleles; ++a) {
            if ((a * nSamples * SIM_PLOIDY) == allelesum) {
                return(-2);
            }
        }
    }

#if 0
    for (int i = 0;i < n_alleles;++i) {
        fprintf(stderr, "Found %s allele: %s with int value %d\n", i == 0 ? "REF" : (i == 1 ? "ALT_1" : (i == 2 ? "ALT_2" : (i == 3 ? "ALT_3" : "ALT_4"))), sim->rec->d.allele[i], rec_alleles[i]);
    }
#endif


    return(0);

}



inline void write_record_values(htsFile* out_fp, simRecord* sim) {

    int ret;

    if (NULL == sim->rec) {
        // last flush at the program end
        ret = prepare_gvcf_block(sim, sim->gvcfd);
        DEVASSERT(GVCF_WRITE_SIMREC != ret);
        if (GVCF_FLUSH_BLOCK == ret) {
            ASSERT(0 == bcf_write(out_fp, sim->hdr, sim->gvcfd->grec));
        }
        return;
    }

    if (NULL == sim->gvcfd) {
        ASSERT(0 == bcf_write(out_fp, sim->hdr, sim->rec));
        return;
    } else {
        // for flushing gvcf blocks (if any)
        // e.g. we have a gvcf block in the buffer but the next record is variant
        // so we need to flush the block then write the variant record
        ret = prepare_gvcf_block(sim, sim->gvcfd);
        if (GVCF_FLUSH_BLOCK == ret) {
            // flush block 
            ASSERT(0 == bcf_write(out_fp, sim->hdr, sim->gvcfd->grec));
            // now try again with current record
            ret = prepare_gvcf_block(sim, sim->gvcfd);
        }
        ASSERT(ret != GVCF_FLUSH_BLOCK); // cannot be asking to flush twice
        if (GVCF_WRITE_SIMREC == ret) {
            ASSERT(0 == bcf_write(out_fp, sim->hdr, sim->rec));
        }
        if (GVCF_NO_WRITE == ret) {
            // write nothing; current block is still growing
            return;
        }
        return;
    }

}


// rng1_seeder: seed for poisson distribution
// used if error_qs != 0 so that error_qs 0 gives the same output as angsd
// isolates poisson distribution RNG from other RNGs
// to be able to sample depths from poisson the same way given the same seed
// regardless of --error-qs beta sampling being on or off
unsigned short int rng1_seeder[3] = SEEDER_INIT;
unsigned short int rng1_seeder_save[3] = SEEDER_INIT;

// rng2_seeder: seed for beta distribution
unsigned short int rng2_seeder[3] = SEEDER_INIT;
unsigned short int rng2_seeder_save[3] = SEEDER_INIT;


void(*calculate_gls)(simRecord* sim);


/// @brief simulate a position with no reads for any of the individuals
/// @param sim 
/// @return 
inline int simulate_site_with_no_reads(simRecord* sim) {

    if (args->printPileup) {
        ksprintf(sim->pileup, "%s\t%ld\t%c", sim->hdr->id[BCF_DT_CTG][sim->rec->rid].key, sim->rec->pos + 1, sim->rec->d.allele[0][0]);
        for (int s = 0; s < sim->nSamples; s++) {
            ksprintf(sim->pileup, "\t0\t*\t*");
        }
        kputc('\n', sim->pileup);
    }

    // remove genotypes from main output file
    bcf_update_genotypes(sim->hdr, sim->rec, NULL, 0);

    // -- gVCF --
    if (args->doGVCF) {
        ASSERT(0 == (bcf_update_alleles_str(sim->hdr, sim->rec, "<NON_REF>")));
        sim->add_tags();
        return (0);
    }

    // -- VCF --

    if (args->doUnobserved == ARG_DOUNOBSERVED_TRIM) {
        ASSERT(0 == (bcf_update_alleles_str(sim->hdr, sim->rec, "."))); //DRAGON missing allele
        sim->nAlleles = 1;
        sim->nGenotypes = 1;
        sim->nAllelesObserved = 0; //DRAGON
    } else if (args->doUnobserved == ARG_DOUNOBSERVED_STAR) {
        ASSERT(0 == (bcf_update_alleles_str(sim->hdr, sim->rec, "<*>")));
        sim->nAlleles = 1;
        sim->nGenotypes = 1;
        sim->nAllelesObserved = 0; //DRAGON
    } else if (args->doUnobserved == ARG_DOUNOBSERVED_NONREF) {
        ASSERT(0 == (bcf_update_alleles_str(sim->hdr, sim->rec, "<NON_REF>")));
        sim->nAlleles = 1;
        sim->nGenotypes = 1;
        sim->nAllelesObserved = 0; //DRAGON
    } else if (args->doUnobserved == ARG_DOUNOBSERVED_EXPLODE_ACGT) {
        ASSERT(0 == (bcf_update_alleles_str(sim->hdr, sim->rec, "A,C,G,T")));
        sim->nAlleles = 4;
        sim->nAllelesObserved = 4;
        sim->nGenotypes = 10;
    } else if (args->doUnobserved == ARG_DOUNOBSERVED_EXPLODE_ACGT_STAR) {
        ASSERT(0 == (bcf_update_alleles_str(sim->hdr, sim->rec, "A,C,G,T,<*>")));
        sim->nAlleles = 5;
        sim->nAllelesObserved = 4;
        sim->nGenotypes = 15;
    } else if (args->doUnobserved == ARG_DOUNOBSERVED_EXPLODE_ACGT_NONREF) {
        ASSERT(0 == (bcf_update_alleles_str(sim->hdr, sim->rec, "A,C,G,T,<NON_REF>")));
        sim->nAlleles = 5;
        sim->nAllelesObserved = 4;
        sim->nGenotypes = 15;
    }

    sim->current_size_bcf_tag_number[FMT_NUMBER_G] = sim->nSamples * sim->nGenotypes;
    sim->current_size_bcf_tag_number[FMT_NUMBER_R] = sim->nSamples * sim->nAllelesObserved;
    sim->current_size_bcf_tag_number[FMT_NUMBER_R_WITH_NONREF] = sim->nSamples * sim->nAlleles;
    sim->current_size_bcf_tag_number[INFO_NUMBER_G] = sim->nGenotypes;
    sim->current_size_bcf_tag_number[INFO_NUMBER_R] = sim->nAllelesObserved;
    sim->current_size_bcf_tag_number[INFO_NUMBER_R_WITH_NONREF] = sim->nAlleles;

    DEVASSERT(sim->current_size_bcf_tag_number[FMT_NUMBER_G] <= sim->max_size_bcf_tag_number[FMT_NUMBER_G]);
    DEVASSERT(sim->current_size_bcf_tag_number[FMT_NUMBER_R] <= sim->max_size_bcf_tag_number[FMT_NUMBER_R]);
    DEVASSERT(sim->current_size_bcf_tag_number[FMT_NUMBER_R_WITH_NONREF] <= sim->max_size_bcf_tag_number[FMT_NUMBER_R_WITH_NONREF]);
    DEVASSERT(sim->current_size_bcf_tag_number[INFO_NUMBER_G] <= sim->max_size_bcf_tag_number[INFO_NUMBER_G]);
    DEVASSERT(sim->current_size_bcf_tag_number[INFO_NUMBER_R] <= sim->max_size_bcf_tag_number[INFO_NUMBER_R]);
    DEVASSERT(sim->current_size_bcf_tag_number[INFO_NUMBER_R_WITH_NONREF] <= sim->max_size_bcf_tag_number[INFO_NUMBER_R_WITH_NONREF]);


    for (int i = 0;i < sim->max_size_bcf_tag_number[bcf_tags[GL].n];++i) {
        if (NULL != sim->pl_arr) {
            sim->pl_arr[i] = bcf_int32_missing;
        }
        if (NULL != sim->gp_arr) {
            sim->gp_arr[i] = bcf_float_missing_union_f;
        }
        sim->gl_arr[i] = bcf_float_missing_union_f;
    }
    // no need to set the rest of the tags to missing since the initialized values (and vals we reset them to during record reset) are 0
    sim->add_tags();

    if (args->printPileup) {
        FLUSH_BGZF_KSTRING_BUFFER(args->out_pileup_fp, sim->pileup);
    }

    return (0);

}






// return value < 0		 skip the site
// returns negative value only if program will skip the position given arg values
// skip position reasons:
// -3	observed only 1 simulated allele at position (iff PROGRAM_WILL_SKIP_SIM_INVAR_SITES)
// -4	no alleles were observed at position (iff 1==args->rmEmptySites)
int simulate_record_values(simRecord* sim) {

    bcf1_t* rec = sim->rec;

    if (0 == rec->n_allele) {
        ERROR("No alleles found at position %ld.", rec->pos + 1);
    }

    if (PROGRAM_WILL_SKIP_INPUT_HOMOGT_SITES) {
        if (1 == rec->n_allele) {
            return(-1);
        }
    }

    const int nSamples = sim->nSamples;
    int n_sim_reads = 0;

    // -- for error-qs 2
    double error_prob_forQs_i = -1.0;
    int qScore_i = -1;
    int adjqScore_i = -1;
    // --

    // -------------------------------------------- //
    // reset reused objects for the current rec
    sim->reset_rec_objects();
    // -------------------------------------------- //


    int s = -1; // sample index
    int b = -1; // base index in ACGT
    int a = -1; // allele index

    // -------------------------------------------- //
    // simulate read depths


    if (args->mps_depths != NULL) {
        poissonSampler_sample_depths_perSample_means(args->poissonSampler, n_sim_reads_arr, nSamples);
    } else {
        poissonSampler_sample_depths_same_mean(args->poissonSampler[0], n_sim_reads_arr, nSamples);
    }

    for (s = 0; s < nSamples; s++) {

        n_sim_reads = n_sim_reads_arr[s];

        sim->fmt_dp_arr[s] = n_sim_reads;
        sim->info_dp_arr[0] += n_sim_reads;

        if (n_sim_reads > sim->_nBasesPerSample) {
            sim->expand_arrays(n_sim_reads);
        }
    }
    // end simulate read depths
    // -------------------------------------------- //

    // -------------------------------------------------------------------- //
    // INFO/DP == 0 //
    // no reads were simulated for any of the individuals //
    if (0 == sim->info_dp_arr[0]) {

        // [FILTER] --rm-empty-sites 
        VWARN("No alleles were observed at position %ld.", sim->rec->pos + 1);
        if (1 == args->rmEmptySites) {
            return (-4);
        }
        return(simulate_site_with_no_reads(sim));
    }

    // -------------------------------------------------------------------- //
    // INFO/DP >0 //

    int which_strand = -1;
    int true_base = -1;  // true base (no error)
    int r_base = -1;     // simulated base (observed base after error)
    int tail_dist = -1;

    if (args->printPileup) {
        ksprintf(sim->pileup, "%s\t%ld\t%c", sim->hdr->id[BCF_DT_CTG][rec->rid].key, rec->pos + 1, rec->d.allele[0][0]);
    }

    int32_t* sample_acgt_fmt_ad_arr = NULL;
    int32_t* sample_acgt_fmt_adf_arr = NULL;
    int32_t* sample_acgt_fmt_adr_arr = NULL;
    int32_t* sample_acgt_fmt_qsum_arr = NULL;
    int32_t* sample_acgt_fmt_qsum_sq_arr = NULL;

    double base_pick_error_prob = args->base_pick_error_prob;
    if (1 == args->error_qs) {
        // if error_qs 1, args->base_pick_error_prob initted to -1.0
        base_pick_error_prob = args->betaSampler->sample();

        if (args->printBasePickError) {
            // TSV: type, sample_id, contig, site, read_index, base_pick_error_prob
            for (s = 0; s < nSamples; s++) {
                fprintf(stdout, "base_pick_error_prob\t%s\t%s\t%ld\tNA\t%f\n", sim->hdr->samples[s], sim->hdr->id[BCF_DT_CTG][rec->rid].key, rec->pos + 1, base_pick_error_prob);
            }
        }

    }


    int* sample_true_gts = NULL;

    for (s = 0; s < nSamples; s++) {

        n_sim_reads = sim->fmt_dp_arr[s];

        if (0 == n_sim_reads) {

            if (args->printPileup) {
                ksprintf(sim->pileup, "\t0\t*\t*");
            }

        } else {

            sample_true_gts = true_gts_acgt_int + (s * SIM_PLOIDY);

            sample_acgt_fmt_ad_arr = sim->acgt_fmt_ad_arr + (s * 4);
            if (1 == args->addFormatADF || 1 == args->addInfoADF) {
                sample_acgt_fmt_adf_arr = sim->acgt_fmt_adf_arr + (s * 4);
            }
            if (1 == args->addFormatADR || 1 == args->addInfoADR) {
                sample_acgt_fmt_adr_arr = sim->acgt_fmt_adr_arr + (s * 4);
            }
            sample_acgt_fmt_qsum_arr = sim->acgt_fmt_qsum_arr + (s * 4);
            sample_acgt_fmt_qsum_sq_arr = sim->acgt_fmt_qsum_sq_arr + (s * 4);

            for (int read_i = 0; read_i < n_sim_reads; read_i++) {

                // -------------------------------------------- //
                // ----> pick a haplotype 
                (sample_uniform_rng1() < 0.5) ? true_base = sample_true_gts[0] : true_base = sample_true_gts[1];


                // -------------------------------------------- //
                // ----> base picking

                // use an indepentent rng for the error base choosing sample_uniform function
                // sample_uniform() in error_base choosing is called for x times where x depends on base_pick_error_prob
                // to be able to sample the same values from haplotype picking and if (sample_uniform() < base_pick_error_prob) statement
                // isolate the base_pick_error_prob dependent part of the RNG from the rest of the RNGs
                // so that we sample the same depths and haplotypes given the same seed when we use error-qs 0, 1 and 2
                // for use with error-qs 1
                r_base = true_base;
                if (sample_uniform_rng0() < base_pick_error_prob) {
                    while ((r_base = (floor(4 * sample_uniform_rng0()))) == true_base);
                }


                // -------------------------------------------- //
                // ----> get qScores

                if (2 == args->error_qs) {
                    error_prob_forQs_i = args->betaSampler->sample();

                    qScore_i = -1;
                    adjqScore_i = -1;

                    if (0.0 == error_prob_forQs_i) {
                        qScore_i = CAP_BASEQ;
                    } else if (1.0 == error_prob_forQs_i) {
                        qScore_i = 0;
                    } else if ((0.0 < error_prob_forQs_i) && (error_prob_forQs_i < 1.0)) {
                        double tmp = -10.0 * log10(error_prob_forQs_i);
                        qScore_i = (int)(tmp);
                        adjqScore_i = (PROGRAM_WILL_ADJUST_QS) ? (int)(tmp + args->adjustBy) : adjqScore_i;
                    } else {
                        ERROR("Bad error probability value: %f", error_prob_forQs_i);
                    }

                    qScore_i = (qScore_i > CAP_BASEQ) ? CAP_BASEQ : qScore_i;
                    qScore_i = (ARG_PLATFORM_RTA3 == args->platform) ? APPLY_RTA3_QSCORE_BINNING(qScore_i) : qScore_i;

                    if (PROGRAM_WILL_ADJUST_QS) {
                        adjqScore_i = (adjqScore_i > CAP_BASEQ) ? CAP_BASEQ : adjqScore_i;
                        adjqScore_i = (ARG_PLATFORM_RTA3 == args->platform) ? (APPLY_RTA3_QSCORE_BINNING(adjqScore_i)) : adjqScore_i;
                    }

                    sim->base_qScores[s][read_i] = qScore_i;

                    if (PROGRAM_WILL_ADJUST_QS) {
                        DEVASSERT(sim->adj_base_qScores != NULL);
                        sim->adj_base_qScores[s][read_i] = adjqScore_i;
                    }


                    if (args->printQsError) {
                        // TSV: type, sample_id, contig, site, read_index, error_prob
                        fprintf(stdout, "qs_error_prob\t%s\t%s\t%ld\t%d\t%f\n", sim->hdr->samples[s], sim->hdr->id[BCF_DT_CTG][rec->rid].key, rec->pos + 1, read_i, error_prob_forQs_i);
                    }

                    if (args->printQScores) {
                        // TSV: type, sample_id, contig, site, read_index, qScore
                        fprintf(stdout, "qs\t%s\t%s\t%ld\t%d\t%d\n", sim->hdr->samples[s], sim->hdr->id[BCF_DT_CTG][rec->rid].key, rec->pos + 1, read_i, ((PROGRAM_WILL_ADJUST_QS_FOR_PRINTQSCORES) ? adjqScore_i : qScore_i));
                    }

                    if (args->usePreciseGlError) {
                        sim->base_error_probs[s][read_i] = error_prob_forQs_i;
                        if (args->printGlError) {
                            // TSV: type, sample_id, contig, site, read_index, error_prob
                            fprintf(stdout, "gl_error_prob\t%s\t%s\t%ld\t%d\t%f\n", sim->hdr->samples[s], sim->hdr->id[BCF_DT_CTG][rec->rid].key, rec->pos + 1, read_i, error_prob_forQs_i);
                        }
                    } else {
                        if (args->printGlError) {
                            // TSV: type, sample_id, contig, site, read_index, error_prob
                            fprintf(stdout, "gl_error_prob\t%s\t%s\t%ld\t%d\t%f\n", sim->hdr->samples[s], sim->hdr->id[BCF_DT_CTG][rec->rid].key, rec->pos + 1, read_i, ((PROGRAM_WILL_ADJUST_QS_FOR_PRINTGLERROR) ? (QS_TO_ERRPROB(adjqScore_i)) : (QS_TO_ERRPROB(qScore_i))));
                        }
                    }


                    if (PROGRAM_WILL_ADJUST_QS_FOR_QSUM) {
                        ASSERT(adjqScore_i != -1);
                        sample_acgt_fmt_qsum_arr[r_base] += adjqScore_i;
                        sample_acgt_fmt_qsum_sq_arr[r_base] += QS_TO_QSSQ(adjqScore_i);
                    } else {
                        sample_acgt_fmt_qsum_arr[r_base] += qScore_i;
                        sample_acgt_fmt_qsum_sq_arr[r_base] += QS_TO_QSSQ(qScore_i);
                    }

                } else {
                    // will use precalculated values at args->preCalc->qScore and args->preCalc->error_prob_forGl instead
                    if (PROGRAM_WILL_ADJUST_QS_FOR_QSUM) {
                        ASSERT(args->preCalc->adj_qScore != -1);
                        sample_acgt_fmt_qsum_arr[r_base] += args->preCalc->adj_qScore;
                        sample_acgt_fmt_qsum_sq_arr[r_base] += QS_TO_QSSQ(args->preCalc->adj_qScore);
                    } else {
                        sample_acgt_fmt_qsum_arr[r_base] += args->preCalc->qScore;
                        sample_acgt_fmt_qsum_sq_arr[r_base] += QS_TO_QSSQ(args->preCalc->qScore);
                    }

                }

                sample_acgt_fmt_ad_arr[r_base]++;

                if (PROGRAM_WILL_SAMPLE_STRAND) {
                    if (sample_uniform_rng0() < 0.5) {
                        which_strand = SIM_FORWARD_STRAND;
                    } else {
                        which_strand = SIM_REVERSE_STRAND;
                    }

                    if (SIM_FORWARD_STRAND == which_strand) {
                        if (sample_acgt_fmt_adf_arr != NULL) {
                            sample_acgt_fmt_adf_arr[r_base]++;
                        }
                    } else if (SIM_REVERSE_STRAND == which_strand) {
                        if (sample_acgt_fmt_adr_arr != NULL) {
                            sample_acgt_fmt_adr_arr[r_base]++;
                        }
                    } else {
                        NEVER;
                    }

                } else {
                    which_strand = SIM_FORWARD_STRAND;
                    if (sample_acgt_fmt_adf_arr != NULL) {
                        sample_acgt_fmt_adf_arr[r_base]++;
                    }
                }

                sim->acgt_n_bases_forI16[(2 * r_base) + which_strand]++;

                sim->bases[s][read_i] = r_base;


            } // read loop


            if (args->printPileup) {
                ksprintf(sim->pileup, "\t%d", n_sim_reads);
                kputc('\t', sim->pileup);
                for (int read_i = 0; read_i < n_sim_reads; read_i++) {
                    kputc("ACGT"[sim->bases[s][read_i]], sim->pileup);
                }
                kputc('\t', sim->pileup);

                if (NULL != sim->base_qScores) {
                    for (int read_i = 0; read_i < n_sim_reads; read_i++) {
                        kputc((PROGRAM_WILL_ADJUST_QS_FOR_PILEUP ? sim->adj_base_qScores[s][read_i] : sim->base_qScores[s][read_i]) + QSCORE_PHRED_ENCODING_OFFSET, sim->pileup);
                    }
                } else {

                    for (int read_i = 0; read_i < n_sim_reads; read_i++) {
                        kputc((PROGRAM_WILL_ADJUST_QS_FOR_PILEUP ? args->preCalc->adj_qScore : args->preCalc->qScore) + QSCORE_PHRED_ENCODING_OFFSET, sim->pileup);
                    }
                }
            }


            for (b = 0; b < 4;++b) {
                sim->acgt_info_ad_arr[b] += sample_acgt_fmt_ad_arr[b];
            }

        }
    } // samples loop

    if (args->printPileup) {
        kputc('\n', sim->pileup);
    }


    if (args->addI16) {

        for (s = 0;s < nSamples;s++) {
            n_sim_reads = sim->fmt_dp_arr[s];
            if (0 != n_sim_reads) {
                for (int read_i = 0; read_i < n_sim_reads; read_i++) {
                    tail_dist = sample_from_range_rng_rand(1, 50);
                    if (tail_dist > CAP_TAIL_DIST) {
                        tail_dist = CAP_TAIL_DIST;
                    }
                    sim->acgt_sum_taildist[r_base] += tail_dist;
                    sim->acgt_sum_taildist_sq[r_base] += (tail_dist * tail_dist);
                }
            }
        }

    }

    int nObservedBases = 0;
    for (b = 0; b < 4; ++b)
    {
        DEVASSERT(sim->acgt_info_ad_arr != NULL);
        if (sim->acgt_info_ad_arr[b] > 0)
        {
            nObservedBases++;
        }
    }

    if (PROGRAM_WILL_SKIP_SIM_INVAR_SITES) {
        if (1 == nObservedBases) {
            return(-3);
        } else if (0 == nObservedBases) {
            NEVER;
        }
    }

    // -------------------------------------------- //
    // -> Reorder alleles: sort by INFO/AD (descending order)
    // use:
    //   acgt_info_ad_arr: contains per-base read depths summed across samples 
    //
    // set: simRecord values:
    //   acgt2alleles
    //   alleles2acgt
    //   nAlleles
    //   nGenotypes
    //   allele_unobserved
    //   alleles

    // step 1) sorting
    // assumption: MAX_NALLELES == 5
    // alelles[4] is reserved for NON_REF
    DEVASSERT(5 == MAX_NALLELES);
    int sorted_indices[4] = { 0,1,2,3 };

    int tmp = -1;
    int n_exploded_bases = 0;

    // insertion sort (descending order)
    for (int i = 1; i < 4; i++) {
        for (int j = i; j > 0 && sim->acgt_info_ad_arr[sorted_indices[j]] > sim->acgt_info_ad_arr[sorted_indices[j - 1]]; j--) {
            tmp = sorted_indices[j];
            sorted_indices[j] = sorted_indices[j - 1];
            sorted_indices[j - 1] = tmp;
        }
    }

    for (a = 0;a < 4;++a) {
        b = sorted_indices[a];
        sim->alleles2acgt[a] = b;
        sim->acgt2alleles[b] = a;
    }

    // step 2) handle unobserved alleles

    int n_observed_bases = 0;
    for (b = 0;b < 4;++b) {
        if (sim->acgt_info_ad_arr[b] > 0) {
            n_observed_bases++;
        } else {
            if (PROGRAM_WILL_EXPLODE_ACGT) {
                // explode A,C,G,T
                ++n_exploded_bases;
            } else {
                sim->alleles2acgt[sim->acgt2alleles[b]] = -1;
                sim->acgt2alleles[b] = -1;
            }
        }
    }

    sim->allele_unobserved = -1;
    int n_alleles = 0;
    for (a = 0; a < 5; ++a) {
        if (-1 == sim->alleles2acgt[a]) {
            if (PROGRAM_WILL_ADD_UNOBSERVED) {
                sim->allele_unobserved = a;
            }
            break;
        }
        ASSERT(a != 4);
        ++n_alleles;
        kputc("ACGT"[sim->alleles2acgt[a]], &sim->alleles);
        if (a != 3 && sim->alleles2acgt[a + 1] != -1) {
            kputc(',', &sim->alleles);
        }
    }


    int n_unobserved = 0;
    if (PROGRAM_WILL_ADD_UNOBSERVED) {
        kputc(',', &sim->alleles);
        kputs(nonref_str, &sim->alleles);
        sim->alleles2acgt[sim->allele_unobserved] = BASE_NONREF;
        sim->acgt2alleles[BASE_NONREF] = sim->allele_unobserved;
        ++n_unobserved;
    }

    sim->nAllelesObserved = n_alleles;
    sim->nAlleles = n_alleles + n_unobserved;
    sim->nGenotypes = NALLELES_TO_NGTS(sim->nAlleles);

    sim->current_size_bcf_tag_number[FMT_NUMBER_G] = sim->nSamples * sim->nGenotypes;
    sim->current_size_bcf_tag_number[FMT_NUMBER_R] = sim->nSamples * sim->nAllelesObserved;
    sim->current_size_bcf_tag_number[FMT_NUMBER_R_WITH_NONREF] = sim->nSamples * sim->nAlleles;
    sim->current_size_bcf_tag_number[INFO_NUMBER_G] = sim->nGenotypes;
    sim->current_size_bcf_tag_number[INFO_NUMBER_R] = sim->nAllelesObserved;
    sim->current_size_bcf_tag_number[INFO_NUMBER_R_WITH_NONREF] = sim->nAlleles;

    DEVASSERT(sim->current_size_bcf_tag_number[FMT_NUMBER_G] <= sim->max_size_bcf_tag_number[FMT_NUMBER_G]);
    DEVASSERT(sim->current_size_bcf_tag_number[FMT_NUMBER_R] <= sim->max_size_bcf_tag_number[FMT_NUMBER_R]);
    DEVASSERT(sim->current_size_bcf_tag_number[FMT_NUMBER_R_WITH_NONREF] <= sim->max_size_bcf_tag_number[FMT_NUMBER_R_WITH_NONREF]);
    DEVASSERT(sim->current_size_bcf_tag_number[INFO_NUMBER_G] <= sim->max_size_bcf_tag_number[INFO_NUMBER_G]);
    DEVASSERT(sim->current_size_bcf_tag_number[INFO_NUMBER_R] <= sim->max_size_bcf_tag_number[INFO_NUMBER_R]);
    DEVASSERT(sim->current_size_bcf_tag_number[INFO_NUMBER_R_WITH_NONREF] <= sim->max_size_bcf_tag_number[INFO_NUMBER_R_WITH_NONREF]);

    ASSERT(0 == (bcf_update_alleles_str(sim->hdr, sim->rec, sim->alleles.s)));

    DEVASSERT(sim->nAlleles == sim->rec->n_allele);

    // -------------------------------------------- //

    calculate_gls(sim);

    // -------------------------------------------- //

    // remove genotypes from main output file
    bcf_update_genotypes(sim->hdr, sim->rec, NULL, 0);


    // -------------------------------------------- //
    // -> set sorted values for *_NUMBER_R_WITH_NONREF tags

    sample_acgt_fmt_ad_arr = NULL;
    sample_acgt_fmt_adf_arr = NULL;
    sample_acgt_fmt_adr_arr = NULL;
    int32_t* sample_fmt_ad_arr = NULL;
    int32_t* sample_fmt_adf_arr = NULL;
    int32_t* sample_fmt_adr_arr = NULL;

    for (s = 0; s < sim->nSamples; ++s) {

        for (a = 0; a < sim->nAlleles; ++a) {

            b = sim->alleles2acgt[a];

            if (BASE_NONREF == b) {
                continue;
            }

            sample_fmt_ad_arr = sim->fmt_ad_arr + (s * sim->nAlleles);
            sample_fmt_adf_arr = sim->fmt_adf_arr + (s * sim->nAlleles);
            sample_fmt_adr_arr = sim->fmt_adr_arr + (s * sim->nAlleles);
            sample_acgt_fmt_ad_arr = sim->acgt_fmt_ad_arr + (s * 4);
            sample_acgt_fmt_adf_arr = sim->acgt_fmt_adf_arr + (s * 4);
            sample_acgt_fmt_adr_arr = sim->acgt_fmt_adr_arr + (s * 4);

            if (NULL != sim->fmt_ad_arr) {
                sample_fmt_ad_arr[a] = sample_acgt_fmt_ad_arr[b];
            }
            if (NULL != sim->fmt_adf_arr) {
                sample_fmt_adf_arr[a] = sample_acgt_fmt_adf_arr[b];
            }
            if (NULL != sim->fmt_adr_arr) {
                sample_fmt_adr_arr[a] = sample_acgt_fmt_adr_arr[b];
            }

            if (NULL != sim->info_ad_arr) {
                sim->info_ad_arr[a] += sample_acgt_fmt_ad_arr[b];
            }
            if (NULL != sim->info_adf_arr) {
                sim->info_adf_arr[a] += sample_acgt_fmt_adf_arr[b];
            }
            if (NULL != sim->info_adr_arr) {
                sim->info_adr_arr[a] += sample_acgt_fmt_adr_arr[b];
            }
        }
    }

    if (1 == args->addQS) {

        // example: at position 0
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
        a = -1;
        b = -1;
        float sum = 0.0;

        sample_acgt_fmt_qsum_arr = NULL;
        sample_acgt_fmt_qsum_sq_arr = NULL;

        // -- calculation of the qsum -- //
        // sum the normalized qsum across all samples
        // to account for differences in coverage
        // modified from: bcftools/bam2bcf.c
        for (s = 0; s < nSamples; ++s) {
            sum = 0.0;
            sample_acgt_fmt_qsum_arr = sim->acgt_fmt_qsum_arr + (s * 4);
            sample_acgt_fmt_qsum_sq_arr = sim->acgt_fmt_qsum_sq_arr + (s * 4);
            for (b = 0; b < 4; ++b) {
                sum += sample_acgt_fmt_qsum_arr[b];
            }

            if (0.0 != sum) {


                for (b = 0; b < 4; ++b) {
                    a = sim->acgt2alleles[b];
                    if (-1 == a) {
                        continue;
                    }

                    // per-sample normalization
                    // save as ordered (base b --> allele index a)
                    sim->qs_arr[a] += (float)((float)(sample_acgt_fmt_qsum_arr[b]) / sum);
                }
            }
        }
    }



    // -------------------------------------------- //
    // -> set sorted values for FMT_NUMBER_G tags
    // 		//GL
    //		GP
    //		PL
    if (args->addPL) {

        int x;
        for (int i = 0; i < sim->current_size_bcf_tag_number[bcf_tags[PL].n];++i) {
            if (bcf_float_is_missing(sim->gl_arr[i])) {
                sim->pl_arr[i] = bcf_int32_missing;
            } else if (bcf_float_is_vector_end(sim->gl_arr[i])) {
                NEVER;
            } else {

#if 0
                if (sim->gl_arr[i] == std::numeric_limits<float>::infinity()) {
                    NEVER;
                }
#endif

                if (sim->gl_arr[i] == NEG_INF) {
                    sim->pl_arr[i] = MAXPL;
                    continue;
                }

                // when max threshold disabled, vcfgl gives same results as bcftools
                // +tag2tag -- --GL-to-PL

                x = lroundf(-10.0 * sim->gl_arr[i]);
                if (x > MAXPL) {
                    x = MAXPL;
                }

                sim->pl_arr[i] = x;
            }
        }
    }

    if (args->addGP) {
        for (int i = 0; i < sim->current_size_bcf_tag_number[bcf_tags[PL].n];++i) {
            if (bcf_float_is_missing(sim->gl_arr[i])) {
                sim->gp_arr[i] = bcf_float_missing_union_f;
            } else if (bcf_float_is_vector_end(sim->gl_arr[i])) {
                NEVER;
            } else {
                sim->gp_arr[i] = pow(10, sim->gl_arr[i]);
            }
        }
    }

    // -------------------------------------------- //
    // prepare I16 tag


// requires:
// acgt_n_bases_forI16
// acgt_fmt_qsum_arr
// acgt_fmt_qsum_sq_arr
// acgt_sum_taildist
// acgt_sum_taildist_sq
    if (1 == args->addI16) {
        // ------------------------------------------------------- //
        // I16 REF fields
        // ------------------------------------------------------- //
        int refb = sim->alleles2acgt[0];

        // 1   #reference Q13 bases on the forward strand
        sim->i16_arr[0] = sim->acgt_n_bases_forI16[(refb * 2) + SIM_FORWARD_STRAND];

        // 2   #reference Q13 bases on the reverse strand
        sim->i16_arr[1] = sim->acgt_n_bases_forI16[(refb * 2) + SIM_REVERSE_STRAND];


        sample_fmt_ad_arr = NULL;

        for (s = 0; s < nSamples; ++s) {
            // 5   sum of reference base qualities
            sim->i16_arr[4] += sim->acgt_fmt_qsum_arr[s * 4 + refb];
            // 6   sum of squares of reference base qualities
            sim->i16_arr[5] += sim->acgt_fmt_qsum_sq_arr[s * 4 + refb];

            sample_fmt_ad_arr = sim->fmt_ad_arr + (s * sim->nAlleles);
            for (a = 0; a < sim->nAlleles; ++a) {
                if (sim->nAllelesObserved == a) {
                    // exclude NONREF allele since that cannot have a mapq
                    continue;
                }
                for (int i = 0; i < sample_fmt_ad_arr[a]; ++i) {
                    if (0 == a) {
                        // 9   sum of ref mapping qualities
                        sim->i16_arr[8] += args->i16_mapq;
                        // 10  sum of squares of ref mapping qualities
                        sim->i16_arr[9] += args->i16_mapq * args->i16_mapq;

                    } else {
                        // 11  sum of non-ref mapping qualities
                        sim->i16_arr[10] += args->i16_mapq;

                        // 12  sum of squares of non-ref mapping qualities
                        sim->i16_arr[11] += args->i16_mapq * args->i16_mapq;
                    }
                }
            }
        }


        // 13  sum of tail distance for ref bases
        sim->i16_arr[12] = sim->acgt_sum_taildist[refb];

        // 14  sum of squares of tail distance for ref bases
        sim->i16_arr[13] = sim->acgt_sum_taildist_sq[refb];

        // ------------------------------------------------------- //
        // I16 NON-REF fields
        // ------------------------------------------------------- //
        // start from 1 to exclude the reference allele
        ASSERT(sim->nAlleles > 1);
        for (a = 1; a < sim->nAlleles; ++a) {

            if (sim->nAllelesObserved == a) {
                continue;
            }

            b = sim->alleles2acgt[a];

            DEVASSERT(b != BASE_NONREF);
            DEVASSERT(b != -1);


            // 3   #non-ref bases on the forward strand
            // 3   (old definition) #non-ref Q13 bases on the forward strand
            sim->i16_arr[2] += sim->acgt_n_bases_forI16[(b * 2) + SIM_FORWARD_STRAND];

            // 4   #non-ref bases on the reverse strand
            // 4   (old definition) #non-ref Q13 bases on the reverse strand
            sim->i16_arr[3] += sim->acgt_n_bases_forI16[(b * 2) + SIM_REVERSE_STRAND];

            for (s = 0; s < nSamples; ++s) {
                // 7   sum of non-ref base qualities
                sim->i16_arr[6] += sim->acgt_fmt_qsum_arr[s * 4 + b];

                // 8   sum of squares of non-ref base qualities
                sim->i16_arr[7] += sim->acgt_fmt_qsum_sq_arr[s * 4 + b];
            }

            // 15  sum of tail distance for non-ref bases
            sim->i16_arr[14] += sim->acgt_sum_taildist[b];

            // 16  sum of squares of tail distance for non-ref
            sim->i16_arr[15] += sim->acgt_sum_taildist_sq[b];
        }

    }


    //

    sim->add_tags();

    if (args->printPileup) {
        FLUSH_BGZF_KSTRING_BUFFER(args->out_pileup_fp, sim->pileup);
    }

    return (0);
}

inline int simulate_record_true_values(simRecord* sim) {

    bcf1_t* rec = sim->rec;
    const int nSamples = sim->nSamples;


    // -------------------------------------------- //
    // reset reused objects for the current rec
    sim->reset_rec_objects();
    // -------------------------------------------- //

    // count of acgt bases in genotypes
    // for A,C,G,T
    int acgt_ac[4] = { 0,0,0,0 };
    for (int i = 0; i < nSamples * SIM_PLOIDY; ++i) {
        acgt_ac[true_gts_acgt_int[i]]++;
    }

    int desc_order_indices[4] = { 0,1,2,3 };
    for (int i = 0;i < 4;++i) {
        if (acgt_ac[i] > 0) {
            sim->nAllelesObserved++;
        }
        for (int j = i;j > 0 && acgt_ac[desc_order_indices[j]] > acgt_ac[desc_order_indices[j - 1]];j--) {
            int tmp = desc_order_indices[j];
            desc_order_indices[j] = desc_order_indices[j - 1];
            desc_order_indices[j - 1] = tmp;
        }
    }

    int k = 0;
    int acgtidx;
    for (int i = 0;i < sim->nAllelesObserved;++i) {
        acgtidx = desc_order_indices[i];
        if (i != 0) {
            kputc(',', &sim->alleles);
        }
        kputc("ACGT"[acgtidx], &sim->alleles);
        ++k;
    }
    if (PROGRAM_WILL_EXPLODE_ACGT) {
        for (int i = sim->nAllelesObserved;i < 4;++i) {
            // desc_order_indices[i]
            acgtidx = desc_order_indices[i];
            kputc(',', &sim->alleles);
            kputc("ACGT"[acgtidx], &sim->alleles);
            ++k;
        }
    }

    if (PROGRAM_WILL_ADD_UNOBSERVED) {
        kputc(',', &sim->alleles);
        kputs(nonref_str, &sim->alleles);
        sim->allele_unobserved = k;
        sim->alleles2acgt[sim->allele_unobserved] = BASE_NONREF;
        sim->acgt2alleles[BASE_NONREF] = sim->allele_unobserved;
    }



    if (PROGRAM_WILL_EXPLODE_ACGT) {

        if (PROGRAM_WILL_ADD_UNOBSERVED) {
            sim->nAlleles = 5;
            sim->nGenotypes = 15;
        } else {
            sim->nAlleles = 4;
            sim->nGenotypes = 10;
        }

    } else {
        // explode 0

        if (PROGRAM_WILL_ADD_UNOBSERVED) {
            sim->nAlleles = sim->nAllelesObserved + 1;
            sim->nGenotypes = NALLELES_TO_NGTS(sim->nAlleles);
        } else {
            sim->nAlleles = sim->nAllelesObserved;
            sim->nGenotypes = NALLELES_TO_NGTS(sim->nAlleles);
        }
    }

    sim->current_size_bcf_tag_number[FMT_NUMBER_G] = sim->nSamples * sim->nGenotypes;
    sim->current_size_bcf_tag_number[FMT_NUMBER_R] = sim->nSamples * sim->nAllelesObserved;
    sim->current_size_bcf_tag_number[FMT_NUMBER_R_WITH_NONREF] = sim->nSamples * sim->nAlleles;
    sim->current_size_bcf_tag_number[INFO_NUMBER_G] = sim->nGenotypes;
    sim->current_size_bcf_tag_number[INFO_NUMBER_R] = sim->nAllelesObserved;
    sim->current_size_bcf_tag_number[INFO_NUMBER_R_WITH_NONREF] = sim->nAlleles;

    DEVASSERT(sim->current_size_bcf_tag_number[FMT_NUMBER_G] <= sim->max_size_bcf_tag_number[FMT_NUMBER_G]);
    DEVASSERT(sim->current_size_bcf_tag_number[FMT_NUMBER_R] <= sim->max_size_bcf_tag_number[FMT_NUMBER_R]);
    DEVASSERT(sim->current_size_bcf_tag_number[FMT_NUMBER_R_WITH_NONREF] <= sim->max_size_bcf_tag_number[FMT_NUMBER_R_WITH_NONREF]);
    DEVASSERT(sim->current_size_bcf_tag_number[INFO_NUMBER_G] <= sim->max_size_bcf_tag_number[INFO_NUMBER_G]);
    DEVASSERT(sim->current_size_bcf_tag_number[INFO_NUMBER_R] <= sim->max_size_bcf_tag_number[INFO_NUMBER_R]);

    ASSERT(0 == (bcf_update_alleles_str(sim->hdr, sim->rec, sim->alleles.s)));
    int true_gt_idx;
    // get 0-based allele indices from the GT tag
    for (int s = 0; s < nSamples; s++) {

        // most likely value should be at the genotype combination for true
        // genotype

        true_gt_idx = bcf_alleles2gt(true_gts_alleles_idx[s * SIM_PLOIDY], true_gts_acgt_int[s * SIM_PLOIDY + 1]);

        for (int i = 0; i < sim->nGenotypes; ++i) {
            if (i == true_gt_idx) {
                sim->gl_arr[s * sim->nGenotypes + i] = MAXGL;

                if (NULL != sim->pl_arr) {
                    sim->pl_arr[s * sim->nGenotypes + i] = MAXPL;
                }
                if (NULL != sim->gp_arr) {
                    sim->gp_arr[s * sim->nGenotypes + i] = MAXGP;
                }
            } else {
                sim->gl_arr[s * sim->nGenotypes + i] = MINGL;

                if (NULL != sim->pl_arr) {
                    sim->pl_arr[s * sim->nGenotypes + i] = MINPL;
                }
                if (NULL != sim->gp_arr) {
                    sim->gp_arr[s * sim->nGenotypes + i] = MINGP;
                }
            }
        }

    }


    // -------------------------------------------- //

    // remove genotypes from main output file
    bcf_update_genotypes(sim->hdr, sim->rec, NULL, 0);
    // -------------------------------------------- //


    if (1 == args->addGL) {
        ASSERT(0 == (bcf_update_format_float(
            sim->hdr, rec, "GL", sim->gl_arr,
            sim->current_size_bcf_tag_number[bcf_tags[GL].n])));
    }

    if (1 == args->addGP) {
        ASSERT(0 == (bcf_update_format_float(
            sim->hdr, rec, "GP", sim->gp_arr,
            sim->current_size_bcf_tag_number[bcf_tags[GP].n])));
    }

    if (1 == args->addPL) {
        ASSERT(0 == (bcf_update_format_int32(
            sim->hdr, rec, "PL", sim->pl_arr,
            sim->current_size_bcf_tag_number[bcf_tags[PL].n])));
    }


    return (0);
}



inline void main_simulate_record_true_values(simRecord* sim, bcf_hdr_t* in_hdr, bcf1_t* in_rec) {

    int nSites = 0;
    int nSitesSkipped = 0;
    int nSitesTotal = 0;
    int ret;

    bcf1_t* explode_rec = NULL;


    while (0 == bcf_read(args->in_fp, in_hdr, in_rec)) {

        bcf_unpack(in_rec, BCF_UN_ALL);

        while (1 == args->explode) {

            if (nSitesTotal == in_rec->pos) {
                break; // out of while(1==args->explode) == out of explode
            }

            if (NULL == explode_rec) {
                // prepare a blank explode record template
                // run only once
                explode_rec = bcf_init();
                explode_rec = bcf_copy(explode_rec, in_rec);
                int32_t* tmp_gt_arr = (int*)malloc(sim->nHaplotypes * sizeof(int));
                for (int h = 0; h < sim->nHaplotypes; ++h) {
                    // assume: input vcf contains all phased gts
                    // therefore the blank rec contains all phased homo gts
                    tmp_gt_arr[h] = BCF_GT_PHASED_0;
                }
                ASSERT(0 == (bcf_update_genotypes(sim->hdr, explode_rec, tmp_gt_arr, sim->nHaplotypes)));
                free(tmp_gt_arr);
                tmp_gt_arr = NULL;

            }

            // -- runs for every unobserved position exploding --
            explode_rec->pos = nSitesTotal;

            // set explode record as the record to use in simulation
            sim->rec = bcf_copy(sim->rec, explode_rec);
            bcf_unpack(sim->rec, BCF_UN_ALL);

            ret = check_rec_alleles(sim);
            if (ret < 0) {
                nSitesSkipped++;
                nSitesTotal++;
                continue;
            } else {
                if (args->printTruth) {
                    ASSERT(0 == bcf_write(args->out_truth_fp, sim->truth_hdr, sim->rec));
                }

                ret = simulate_record_true_values(sim);
                if (ret < 0) {
                    nSitesSkipped++;
                    nSitesTotal++;
                    continue;
                }
            }

            ASSERT(0 == bcf_write(args->out_fp, sim->hdr, sim->rec));
            nSites++;
            nSitesTotal++;

        }

        sim->rec = bcf_copy(sim->rec, in_rec);
        bcf_unpack(sim->rec, BCF_UN_ALL);

        ret = check_rec_alleles(sim);
        if (ret < 0) {
            nSitesSkipped++;
            nSitesTotal++;
            continue;
        } else {
            if (args->printTruth) {
                ASSERT(0 == bcf_write(args->out_truth_fp, sim->truth_hdr, sim->rec));
            }

            ret = simulate_record_true_values(sim);
            if (ret < 0) {
                nSitesSkipped++;
                nSitesTotal++;
                continue;
            }
        }



        ASSERT(0 == bcf_write(args->out_fp, sim->hdr, sim->rec));
        nSites++;
        nSitesTotal++;

    }


    // read all the records in the vcf file
    // if explode and not end of contig, then simulate records until the end of the contig

    const int contigsize = in_hdr->id[BCF_DT_CTG][in_rec->rid].val->info[0];
    while (1 == args->explode) {

        if (nSitesTotal == contigsize) {
            bcf_destroy(explode_rec);
            break;
        }

        explode_rec->pos = nSitesTotal;

        sim->rec = bcf_copy(sim->rec, explode_rec);
        bcf_unpack(sim->rec, BCF_UN_ALL);


        ret = check_rec_alleles(sim);
        if (ret < 0) {
            nSitesSkipped++;
            nSitesTotal++;
            continue;
        } else {
            if (args->printTruth) {
                ASSERT(0 == bcf_write(args->out_truth_fp, sim->truth_hdr, sim->rec));
            }

            ret = simulate_record_true_values(sim);
            if (ret < 0) {
                nSitesSkipped++;
                nSitesTotal++;
                continue;
            }
        }




        ASSERT(0 == bcf_write(args->out_fp, sim->hdr, sim->rec));
        nSites++;
        nSitesTotal++;
    }


    // last write for flushing the gvcf block if any
    bcf_destroy(sim->rec); // first destroy the old rec
    sim->rec = NULL;


    // /END/ main sites loop -----------------------------------------------

    fprintf(stderr, "\n\n-> Simulation finished successfully.\n\nSummary:\n\tNumber of samples: %i\n\tNumber of contigs: %d\n\tTotal number of sites simulated: %i\n\tNumber of sites included in simulation output file: %i\n\tNumber of sites skipped: %i\n", sim->nSamples, sim->nContigs, nSitesTotal, nSites, nSitesSkipped);
    fprintf(args->arg_fp, "\n\n-> Simulation finished successfully.\n\nSummary:\n\tNumber of samples: %i\n\tNumber of contigs: %d\n\tTotal number of sites simulated: %i\n\tNumber of sites included in simulation output file: %i\n\tNumber of sites skipped: %i\n", sim->nSamples, sim->nContigs, nSitesTotal, nSites, nSitesSkipped);

    // /END/ main sites loop -----------------------------------------------


}

inline void main_simulate_record_values(simRecord* sim, bcf_hdr_t* in_hdr, bcf1_t* in_rec) {

    int nSites = 0;
    int nSitesSkipped = 0;
    int nSitesTotal = 0;
    int ret;

    bcf1_t* explode_rec = NULL;


    while (0 == bcf_read(args->in_fp, in_hdr, in_rec)) {

        bcf_unpack(in_rec, BCF_UN_ALL);


        while (1 == args->explode) {

            if (nSitesTotal == in_rec->pos) {
                break; // out of while(1==args->explode) == out of explode
            }

            if (NULL == explode_rec) {
                // prepare a blank explode record template
                // run only once
                explode_rec = bcf_init();
                explode_rec = bcf_copy(explode_rec, in_rec);
                int32_t* tmp_gt_arr = (int*)malloc(sim->nHaplotypes * sizeof(int));
                for (int h = 0; h < sim->nHaplotypes; ++h) {
                    // assume: input vcf contains all phased gts
                    // therefore the blank rec contains all phased homo gts
                    tmp_gt_arr[h] = BCF_GT_PHASED_0;
                }
                ASSERT(0 == (bcf_update_genotypes(sim->hdr, explode_rec, tmp_gt_arr, sim->nHaplotypes)));
                free(tmp_gt_arr);
                tmp_gt_arr = NULL;
            }

            // -- runs for every unobserved position exploding --
            explode_rec->pos = nSitesTotal;

            // set explode record as the record to use in simulation
            sim->rec = bcf_copy(sim->rec, explode_rec);
            bcf_unpack(sim->rec, BCF_UN_ALL);


            ret = check_rec_alleles(sim);
            if (ret < 0) {
                nSitesSkipped++;
                nSitesTotal++;
                continue;
            } else {
                if (args->printTruth) {
                    ASSERT(0 == bcf_write(args->out_truth_fp, sim->truth_hdr, sim->rec));
                }

                ret = simulate_record_values(sim);
                if (ret < 0) {
                    nSitesSkipped++;
                    nSitesTotal++;
                    continue;
                }
            }



            write_record_values(args->out_fp, sim);
            nSites++;
            nSitesTotal++;

        }

        sim->rec = bcf_copy(sim->rec, in_rec);
        bcf_unpack(sim->rec, BCF_UN_ALL);


        ret = check_rec_alleles(sim);
        if (ret < 0) {
            nSitesSkipped++;
            nSitesTotal++;
            continue;
        } else {
            if (args->printTruth) {
                ASSERT(0 == bcf_write(args->out_truth_fp, sim->truth_hdr, sim->rec));
            }

            ret = simulate_record_values(sim);
            if (ret < 0) {
                nSitesSkipped++;
                nSitesTotal++;
                continue;
            }
        }



        write_record_values(args->out_fp, sim);
        nSites++;
        nSitesTotal++;

    }


    // we've read all the records in the vcf file
    // if explode and not end of contig, then simulate records until the end of the contig

    const int contigsize = in_hdr->id[BCF_DT_CTG][in_rec->rid].val->info[0];
    while (1 == args->explode) {

        if (nSitesTotal == contigsize) {
            bcf_destroy(explode_rec);
            break;
        }

        explode_rec->pos = nSitesTotal;

        sim->rec = bcf_copy(sim->rec, explode_rec);
        bcf_unpack(sim->rec, BCF_UN_ALL);

        ret = check_rec_alleles(sim);
        if (ret < 0) {
            nSitesSkipped++;
            nSitesTotal++;
            continue;
        } else {
            if (args->printTruth) {
                ASSERT(0 == bcf_write(args->out_truth_fp, sim->truth_hdr, sim->rec));
            }

            ret = simulate_record_values(sim);
            if (ret < 0) {
                nSitesSkipped++;
                nSitesTotal++;
                continue;
            }
        }

        write_record_values(args->out_fp, sim);
        nSites++;
        nSitesTotal++;
    }


    // last write for flushing the gvcf block if any
    bcf_destroy(sim->rec); // first destroy the old rec
    sim->rec = NULL;
    if (NULL != sim->gvcfd) {
        write_record_values(args->out_fp, sim);
    }

    // /END/ main sites loop -----------------------------------------------

    fprintf(stderr, "\n\n-> Simulation finished successfully.\n\nSummary:\n\tNumber of samples: %i\n\tNumber of contigs: %d\n\tTotal number of sites simulated: %i\n\tNumber of sites included in simulation output file: %i\n\tNumber of sites skipped: %i\n", sim->nSamples, sim->nContigs, nSitesTotal, nSites, nSitesSkipped);
    fprintf(args->arg_fp, "\n\n-> Simulation finished successfully.\n\nSummary:\n\tNumber of samples: %i\n\tNumber of contigs: %d\n\tTotal number of sites simulated: %i\n\tNumber of sites included in simulation output file: %i\n\tNumber of sites skipped: %i\n", sim->nSamples, sim->nContigs, nSitesTotal, nSites, nSitesSkipped);

    // /END/ main sites loop -----------------------------------------------


}


int main(int argc, char** argv) {

    clock_t tic = clock();

    args = args_get(--argc, ++argv);

    time_t start_time = time(NULL);
    struct tm* local_start_time = NULL;
    ASSERT((local_start_time = localtime(&start_time)) != NULL);
    ASSERT(strftime(args->datetime, 256, "%a %b %d %H:%M:%S %Y", local_start_time) != 0);

    fprintf(stderr, "\n%s\n\n", args->command);
    fprintf(stderr, "\n[Program start] %s\n", args->datetime);

    fprintf(args->arg_fp, "\n%s\n\n", args->command);
    fprintf(args->arg_fp, "\n[Program start] %s\n", args->datetime);


    if ((0 == args->error_qs) || (1 == args->error_qs)) {
        args->preCalc = new preCalcStruct();
        args->preCalc->error_prob_forQs = args->error_rate;

        int qs = -1;
        int adjqs = -1;

        if (0.0 == args->preCalc->error_prob_forQs) {
            qs = CAP_BASEQ;
            adjqs = CAP_BASEQ;
        } else if (1.0 == args->preCalc->error_prob_forQs) {
            qs = 0;
            adjqs = 0;
        } else if ((0.0 < args->preCalc->error_prob_forQs) && (args->preCalc->error_prob_forQs < 1.0)) {
            double tmp = -10.0 * log10(args->preCalc->error_prob_forQs);
            qs = (int)(tmp);
            if (PROGRAM_WILL_ADJUST_QS) {
                adjqs = (int)(tmp + args->adjustBy);
            }
        } else {
            ERROR("Bad error probability value: %f", args->preCalc->error_prob_forQs);
        }

        qs = (qs > CAP_BASEQ) ? CAP_BASEQ : qs;
        qs = (ARG_PLATFORM_RTA3 == args->platform) ? APPLY_RTA3_QSCORE_BINNING(qs) : qs;
        args->preCalc->qScore = qs;
        args->preCalc->q5 = args->preCalc->qScore << 5;

        if (PROGRAM_WILL_ADJUST_QS) {
            adjqs = (adjqs > CAP_BASEQ) ? CAP_BASEQ : adjqs;
            adjqs = (ARG_PLATFORM_RTA3 == args->platform) ? APPLY_RTA3_QSCORE_BINNING(adjqs) : adjqs;
            args->preCalc->adj_qScore = adjqs;
            args->preCalc->adj_q5 = args->preCalc->adj_qScore << 5;
        }


        if (1 == args->GL) {

            calculate_gls = alleles_calculate_gls_log10_glModel1_fixedQScore;

        } else if (2 == args->GL) {

            if (0 == args->usePreciseGlError) {
                if (PROGRAM_WILL_ADJUST_QS_FOR_GL) {
                    args->preCalc->error_prob_forGl = QS_TO_ERRPROB(args->preCalc->adj_qScore);
                    args->preCalc->homT = qScore_to_log10_gl[0][args->preCalc->adj_qScore];
                    args->preCalc->het = qScore_to_log10_gl[1][args->preCalc->adj_qScore];
                    args->preCalc->homF = qScore_to_log10_gl[2][args->preCalc->adj_qScore];
                } else {
                    args->preCalc->error_prob_forGl = QS_TO_ERRPROB(args->preCalc->qScore);
                    args->preCalc->homT = qScore_to_log10_gl[0][args->preCalc->qScore];
                    args->preCalc->het = qScore_to_log10_gl[1][args->preCalc->qScore];
                    args->preCalc->homF = qScore_to_log10_gl[2][args->preCalc->qScore];
                }

            } else if (1 == args->usePreciseGlError) {

                args->preCalc->error_prob_forGl = args->preCalc->error_prob_forQs;
                if (0.0 == args->preCalc->error_prob_forGl) {
                    args->preCalc->homT = 0;
                    args->preCalc->het = -0.30103;
                    args->preCalc->homF = NEG_INF;
                } else if (args->preCalc->error_prob_forGl > 0.0) {
                    // prepare for gl calculation with pre-calculated terms
                    args->preCalc->homT = log10(1.0 - args->preCalc->error_prob_forGl);
                    args->preCalc->het = log10((1.0 - args->preCalc->error_prob_forGl) / 2.0 + args->preCalc->error_prob_forGl / 6.0);
                    args->preCalc->homF = log10(args->preCalc->error_prob_forGl) - PRE_CALC_LOG10_3;
                } else {
                    NEVER;
                }
            }

            calculate_gls = alleles_calculate_gls_log10_glModel2_fixedQScore;
        }

        if (args->printGlError) {
            fprintf(stdout, "gl_error_prob\tNA\tNA\tNA\tNA\t%f\n", args->preCalc->error_prob_forGl);
        }

        if (args->printQsError) {
            fprintf(stdout, "qs_error_prob\tNA\tNA\tNA\tNA\t%f\n", args->preCalc->error_prob_forQs);
        }

        if (args->printQScores) {
            fprintf(stdout, "qs\tNA\tNA\tNA\tNA\t%d\n", (PROGRAM_WILL_ADJUST_QS) ? args->preCalc->adj_qScore : args->preCalc->qScore);
        }

    } else if (2 == args->error_qs) {
        if (1 == args->GL) {
            calculate_gls = alleles_calculate_gls_log10_glModel1;
        } else if (2 == args->GL) {
            if (args->usePreciseGlError) {
                calculate_gls = alleles_calculate_gls_log10_glModel2_precise1;
            } else {
                calculate_gls = alleles_calculate_gls_log10_glModel2_precise0;
            }
        }
    }


    args->in_fp = open_htsFile(args->in_fn, "r");

    bcf_hdr_t* in_hdr = bcf_hdr_read(args->in_fp);
    bcf1_t* in_rec = bcf_init();
    simRecord* sim = new simRecord(in_hdr);

    true_gts_acgt_int = (int*)malloc(sim->nSamples * SIM_PLOIDY * sizeof(int));
    ASSERT(NULL != true_gts_acgt_int);

    true_gts_alleles_idx = (int*)malloc(sim->nSamples * SIM_PLOIDY * sizeof(int));
    ASSERT(NULL != true_gts_alleles_idx);

    n_sim_reads_arr = (int*)malloc(sim->nSamples * sizeof(int));
    ASSERT(NULL != n_sim_reads_arr);

    args->out_fp = open_htsFile(args->out_fn, args->output_mode_str);
    if (NULL != args->out_truth_fn) {
        args->out_truth_fp = open_htsFile(args->out_truth_fn, args->output_mode_str);
    }
    if (args->printPileup) {
        args->out_pileup_fp = open_BGZF(args->out_pileup_fn, "w");
    }

    // create multithreaded pool
    htsThreadPool tpool = { NULL, 0 };
    if (args->n_threads > 1) {
        tpool.pool = hts_tpool_init(args->n_threads);
        ASSERT(NULL != tpool.pool);
        // add input stream to the pool
        hts_set_opt(args->in_fp, HTS_OPT_THREAD_POOL, &tpool);
        // add output stream to the pool
        hts_set_opt(args->out_fp, HTS_OPT_THREAD_POOL, &tpool);
        if (NULL != args->out_truth_fn) {
            // add truth output stream to the pool
            hts_set_opt(args->out_truth_fp, HTS_OPT_THREAD_POOL, &tpool);
        }
    }

    ASSERT(bcf_hdr_write(args->out_fp, sim->hdr) == 0);
    if (NULL != args->out_truth_fn) {
        ASSERT(bcf_hdr_write(args->out_truth_fp, sim->truth_hdr) == 0);
    }

    if (ARG_DEPTH_INF == args->mps_depth) {
        main_simulate_record_true_values(sim, in_hdr, in_rec);
    } else {
        main_simulate_record_values(sim, in_hdr, in_rec);
    }

    bcf_hdr_destroy(in_hdr);
    bcf_destroy(in_rec);

    ASSERT(0 == hts_close(args->in_fp));
    ASSERT(0 == hts_close(args->out_fp));
    if (NULL != args->out_truth_fp) {
        ASSERT(0 == hts_close(args->out_truth_fp));
    }


    if (args->printPileup) {
        if (sim->pileup->l > 0) {
            write_BGZF(args->out_pileup_fp, sim->pileup->s, sim->pileup->l);
            sim->pileup->l = 0;
        }
        if (NULL != args->out_pileup_fp) {
            ASSERT(0 == bgzf_close(args->out_pileup_fp));
            args->out_pileup_fp = NULL;
        }
    }

    if (NULL != tpool.pool) {
        hts_tpool_destroy(tpool.pool);
    }

    delete sim;

    free(true_gts_acgt_int);
    true_gts_acgt_int = NULL;

    free(true_gts_alleles_idx);
    true_gts_alleles_idx = NULL;

    free(n_sim_reads_arr);
    n_sim_reads_arr = NULL;

    clock_t toc = clock();

    fprintf(stderr, "\n\tElapsed time (CPU): %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
    fprintf(args->arg_fp, "\n\tElapsed time (CPU): %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

    time_t end_time = time(NULL);
    fprintf(stdout, "\tElapsed time (Real): %f seconds\n", difftime(end_time, start_time));
    fprintf(args->arg_fp, "\tElapsed time (Real): %f seconds\n", difftime(end_time, start_time));



    args_destroy(args);

    return(0);
}


