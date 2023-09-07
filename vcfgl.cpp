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

#include "random_generator.h"
#include "estimator.h"
#include "io.h"
#include "version.h"
#include "shared.h"
#include "lut.h"
#include "bcf_utils.h"


/// @brief sample_strand - sample a strand
/// @return		0: forward, 1: reverse
inline int sample_strand(){
	return ( drand48() < 0.5 ? 0 : 1);
}


/// @brief sample_tail_distance - sample a tail length for I16 tag 
/// @return		int length
inline int sample_tail_distance(){
	int i;
	if((i=sample_uniform_from_range(1,50))>CAP_DIST){
		return(CAP_DIST);
	}
	return(i);
}


FILE *getFILE(const char*fname,const char* mode){
	FILE *fp;
	if(NULL==(fp=fopen(fname,mode))){
		fprintf(stderr,"[%s:%s()]\t->Error opening FILE handle for file:%s exiting\n",__FILE__,__FUNCTION__,fname);
		exit(0);
	}
	return fp;
}



FILE *openFile(const char* a,const char* b){
	char *c = (char*)malloc(strlen(a)+strlen(b)+1);
	strcpy(c,a);
	strcat(c,b);
	// fprintf(stderr,"\t-> Dumping file: %s\n",c);
	FILE *fp = getFILE(c,"w");
	free(c);
	return fp;
}



int pick_base(double errate, int inbase){
	int outbase;
	if (drand48()<errate){
		while ((outbase=(floor(4*drand48()))) == inbase);
		return outbase;
	}
	else return inbase;
}


char *get_time(){
	time_t current_time;
	struct tm *local_time; 
	current_time=time(NULL);
	local_time=localtime(&current_time);
	return(asctime(local_time));
}

int32_t *gt_arr=NULL;


int setblank(bcf1_t *blk,bcf1_t *unmod,bcf_hdr_t *hdr){

	int32_t ngt_arr=0;
	int ngt=bcf_get_genotypes(hdr, blk, &gt_arr, &ngt_arr);
	if ( ngt<=0 ){
		fprintf(stderr,"\nGT not present\n");
	}

	int32_t *tmpia = gt_arr;
	for(int i=0; i<bcf_hdr_nsamples(hdr);i++){
		tmpia[2*i+0] = bcf_gt_phased(0);
		tmpia[2*i+1] = bcf_gt_phased(0);
	}
	ASSERT(0==(bcf_update_genotypes(hdr, blk, tmpia, bcf_hdr_nsamples(hdr)*2)));

	return 0;
}






int simulate_record_values(bcf_hdr_t *out_hdr,bcf1_t *out_bcf_rec,int nSamples,double* mps_depths, argStruct* args, const int site_i, FILE* out_baseCounts_ff){

	kstring_t* alleles_str=kbuf_init();

	double* gls= (double*) malloc(nSamples*SIM_NGTS* sizeof(double));
	ASSERT(NULL!=gls);
	for (int i=0;i<nSamples*SIM_NGTS; ++i){
		gls[i]=-0.0;
	}

	int nGenotypes=SIM_NGTS;


	int32_t *dp_vals=NULL;
	dp_vals= bcf_tag_alloc<int32_t>(DP, 0);


	int n_sim_reads=-1; 

	// qScore	phred-scaled quality score
	// 			qScore = -10 * log10(error_probability)
	int qScore=-1;

	int which_strand, which_haplo=-1;
	int refnref=-1;
	int e_base, base=-1;
	double e=-1.0; // simulated error rate instance


	int32_t ngt_arr=0;
	int ngt=bcf_get_genotypes(out_hdr, out_bcf_rec, &gt_arr, &ngt_arr);
	if ( ngt<=0 ){
		ERROR("Could not find GT tag.");
	}

	if(2!=ngt/nSamples){
		ERROR("Ploidy %d is not supported.\n",ngt/nSamples);
	}

	// refnref_n_q13_bases	number of bases where quality score >= 13 summed across all individuals
	//
	// \def refnref_n_q13_bases[2][2]
	// 		refnref_n_q13_bases[reference|non-reference][forward-strand|reverse-strand]
	//
	// e.g.
	// 		refnref_n_q13_bases[0][0] = #reference bases on forward-strand with q >= 13
	int refnref_n_q13_bases[2][2]={{0,0},{0,0}};


	// refnref_sum_taildist	sum of tail lengths
	// e.g.
	// 		refnref_sum_taildist[1]	= sum of tail lengths for non-ref bases
	//
	int refnref_sum_taildist[2]={0,0};
	int refnref_sum_taildist_sq[2]={0,0};


	int refnref_sum_qs[2]={0,0};
	int refnref_sum_qs_sq[2]={0,0};

	// \def acgt_counts[4]	base counts
	int acgt_counts[4]={0};


	// per base qscore sums for each sample
	// 	where alleles = {A, C, G, T} (thus [4])
	int samples_acgt_qs[4][nSamples];
	for(int j=0;j<4;++j){
		for(int i=0;i<nSamples;++i){
			samples_acgt_qs[j][i]=0.0;
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
	int acgt2alleles[5];// alleles: ref, alt, alt2, alt3...
	acgt2alleles[0]=acgt2alleles[1]=acgt2alleles[2]=acgt2alleles[3]=acgt2alleles[4]=-1;


	int tail_dist=0;

	for (int sample_i=0; sample_i<nSamples; sample_i++) {

		if(args->printBaseCounts==1){
			acgt_counts[0]=acgt_counts[1]=acgt_counts[2]=acgt_counts[3]=0;
		}

		if(mps_depths!=NULL){
			n_sim_reads=Poisson(mps_depths[sample_i]);
		}else{
			n_sim_reads=Poisson(args->mps_depth);
		}

		if(n_sim_reads==0){

			dp_vals[sample_i]=0;

		}else{
			

			int32_t *gt_ptr = gt_arr + sample_i*SIM_PLOIDY;

			double* sample_gls=gls + sample_i * SIM_NGTS;

			int bin_gts[2] = {0};


			//binary input genotypes from simulated input
			for (int i=0; i<SIM_PLOIDY;i++){
				bin_gts[i]= bcf_gt_allele(gt_ptr[i]);
				// use bit shifting to check if bin_gt[n] is 0 or 1
				if( (bin_gts[i] >> 1) != 0){
					fprintf(stderr,"ERROR:\n\nbin_gts[%d]: Genotype %d not supported.\n",i,bin_gts[i]);
					exit(1);
				}
			}

			int sim_strands[n_sim_reads]={-1};
			// int sim_quals[n_sim_reads]={-1};


			for (int i=0; i<n_sim_reads; i++){

				(drand48()<0.5) ? which_haplo=0 : which_haplo=1;

				e=drand48();

				e_base = base= bin_gts[which_haplo];
				if(e < args->errate){
					do{
						e_base=floor(4*drand48());
					}while (e_base == base);
					base=e_base;
				}

				// in expected input file, ref is always 0 and non-ref is always 1
				// refnref == 0 -> we have ref
				// refnref == 1 -> we have nonref 
				(0==base) ? refnref=0: refnref=1;

				gl_log10(base,args->errate,sample_gls);


				acgt_counts[base]++;

				qScore=floor(-10 * log10(e));

				(1==args->addI16)? which_strand=sample_strand() : which_strand=1;
				sim_strands[i]=which_strand;

				if(qScore>=13){
					refnref_n_q13_bases[refnref][which_strand]++;
				}

				refnref_sum_qs[refnref]=qScore;
				refnref_sum_qs_sq[refnref]=lut_qscore2[qScore];


				tail_dist=sample_tail_distance();
				refnref_sum_taildist[refnref]+=tail_dist;
				refnref_sum_taildist_sq[refnref]+=tail_dist*tail_dist;

				samples_acgt_qs[base][sample_i]+=qScore;

			}

			rescale_likelihood_ratio(sample_gls);
			dp_vals[sample_i]=n_sim_reads;

		}


		if(NULL!=out_baseCounts_ff){
			fprintf(out_baseCounts_ff,"%d\t%d\t%d\t%d\t%d\t%d\n",site_i,sample_i,acgt_counts[0],acgt_counts[1],acgt_counts[2],acgt_counts[3]);
		}


	} //end sample loop



	
	// [BEGIN QS TAG] ------------------------------------------------------- //
	float* qs_vals=NULL;
	if(1==args->addQS){
	
		qs_vals= bcf_tag_alloc<float>(QS,0.0);

		// format after seq_nt16_str:
		// 0,1,2,3 for A,C,G,T; 4 otherwise. <0 indel
		// int ori_ref=seq_nt16_int[seq_nt16_table[SIM_acgt2alleles[0]]];
		int ori_ref=0; // equivalent to above

		int ref4=0; // TODO 

		int unseen=-1;
		acgt2alleles[0]=ref4;
		int n_alleles=0;

		float qsum[5]={0.0};

		// -- calculation of the qsum -- //
		// sum the normalized qsum across all samples
		// to account for differences in coverage
		// modified from: bcftools/bam2bcf.c
		for(int i=0;i<nSamples;++i){
			
			float sum=0.0;
			for(int j=0;j<4;++j){
				sum += samples_acgt_qs[j][i];
			}

			if(0!=sum){
				for(int j=0;j<4;++j){
					// per-sample normalization
					qsum[j] += (float) samples_acgt_qs[j][i]/sum;
				}
			}else{
				// WARNING("Sum of the quality scores at site with index %d is %d.",site_i);
			}
		}


		// sort qsum in ascending order (insertion sort)
		float *ptr[5], *tmp;

		for (int i=0; i<5; i++){
			ptr[i] = &qsum[i];
		}

		for (int i=1; i<4; i++){
			for (int j=i; j>0 && *ptr[j] < *ptr[j-1]; j--){
				tmp = ptr[j], ptr[j] = ptr[j-1], ptr[j-1] = tmp;
			}
		}

		int i,j=0;
		// Set the reference allele and alternative allele(s)
		for (i=3, j=1; i>=0; i--)   // i: alleles sorted by QS; j, acgt2alleles[j]: output allele ordering
		{
			int ipos = ptr[i] - qsum;   // position in sorted qsum array
			ASSERT(ipos>-1 && ipos<5);
			if ( ipos==ref4 ){
				qs_vals[0] = qsum[ipos];    // REF's qsum
			}
			else
			{
				// NB this will cause the i to be >0 in the below checks
				if ( !qsum[ipos] ) {
					// NEVER;
					break;       // qsum is 0, this and consequent alleles are not seen in the pileup
				}
				qs_vals[j] = qsum[ipos];
				acgt2alleles[j++]  = ipos;
			}
		}
		int ref_base=0; // TODO delme unnec var & ifs nextline
		if (ref_base >= 0)
		{

			// for SNPs, find the "unseen" base
			if (((ref4 < 4 && j < 4) || (ref4 == 4 && j < 5)) && i >= 0){
				// NEVER;
				unseen = j, acgt2alleles[j++] = ptr[i] - qsum;
			}
			n_alleles = j;
		}
		else
		{
			NEVER;
			// n_alleles = j;
			// // if (n_alleles == 1) return -1; // no reliable supporting read. stop doing anything
			// if (n_alleles == 1) NEVER;
		}
		// int has_alt = (n_alleles==2 && unseen!=-1) ? 0 : 1;
		// ASSERT(unseen==-1);
		// ASSERT(has_alt==1);

		nGenotypes=nAlleles2nGenotypes(n_alleles);

		int nals=1;
		kputc("ACGTN"[ori_ref], alleles_str);
		for (int i=1; i<5; ++i)
		{
			if (acgt2alleles[i] < 0) break;
			kputc(',', alleles_str);

			// if ( unseen==i ){
				// // NEVER;
				// kputs("<*>", alleles_str);
			// }else{
				// kputc("ACGT"[acgt2alleles[i]], alleles_str);
				// has_alt = 1;
			// }

			//@@
			kputc("ACGT"[acgt2alleles[i]], alleles_str);
			//@@
				// has_alt = 1;


			nals++;
		}

		
	}

	// [END QS TAG] --------------------------------------------------------- //
	


	float *gl_vals=NULL;
	gl_vals= bcf_tag_alloc<float>(GL,-0.0,nSamples*nGenotypes);


	// set REF and ALT alleles
	if(0==alleles_str->l){
		ASSERT(1!=args->addQS);
		kputs("A,C,G,T", alleles_str);
	}
	// ASSERT(alleles_str->l>1);
	ASSERT(0==(bcf_update_alleles_str(out_hdr,out_bcf_rec,alleles_str->s)));



	ASSERT(0==(bcf_update_format_int32(out_hdr, out_bcf_rec, "DP", dp_vals,bcf_tags[DP].n)));


	// [BEGIN GL TAG] ------------------------------------------------------- //
	for (int sample_i=0;sample_i<nSamples; ++sample_i){

			double* sample_gls=gls + sample_i * SIM_NGTS;

		if(0==dp_vals[sample_i]){
			if(0==args->addQS){
				ASSERT(nGenotypes==SIM_NGTS);
			}
			for(int j=0;j<nGenotypes;j++){
				bcf_float_set_missing(gl_vals[sample_i*nGenotypes+j]);
			}
		}else{
			if(0==args->addQS){
				ASSERT(nGenotypes==SIM_NGTS);
				for(int j=0;j<nGenotypes;j++){
					gl_vals[sample_i*nGenotypes+j]= (float) sample_gls[acgt_genotype_order_lut[j]];
				}
			}else{

				for(int a1=0;a1<4;++a1){
					for(int a2=a1;a2<4;++a2){

						int aa1=acgt2alleles[a1];
						int aa2=acgt2alleles[a2];
						if(aa1==-1 || aa2==-1) continue;

						int newoffset=bcf_alleles2gt(a1,a2);
						int oldoffset=bcf_alleles2gt(aa1,aa2);

						ASSERT(oldoffset>-1);
						ASSERT(newoffset>-1);
						gl_vals[sample_i*nGenotypes+newoffset]= (float) sample_gls[acgt_genotype_order_lut[oldoffset]];

					}
				}
			}
		}

	}

	kbuf_destroy(alleles_str);


	ASSERT(0==(bcf_update_format_float(out_hdr, out_bcf_rec, "GL", gl_vals,bcf_tags[GL].n)));
	// [END GL TAG] --------------------------------------------------------- //


	int32_t *pl_vals=NULL;
	if(1==args->addPL){
		pl_vals = bcf_tag_alloc<int32_t>(PL,0, nSamples*nGenotypes);

		int x;
		for (int i=0; i<bcf_tags[PL].n; ++i)
		{
			if ( bcf_float_is_missing(gl_vals[i]) ){
				pl_vals[i] = bcf_int32_missing;
			} else if ( bcf_float_is_vector_end(gl_vals[i]) ){
				NEVER;
				// pl_vals[i] = bcf_int32_vector_end;
			}
			else {
				((x = lroundf(-10*gl_vals[i])) > MAXPL) ? pl_vals[i]=MAXPL : pl_vals[i]=x;
			}
		}
		ASSERT(0==(bcf_update_format_int32(out_hdr,out_bcf_rec,"PL",pl_vals,bcf_tags[PL].n)));

		free(pl_vals);
		pl_vals=NULL;
	}

	float* gp_vals=NULL;

	if(1==args->addGP){
		gp_vals=bcf_tag_alloc<float>(GP,0.0, nSamples*nGenotypes);
		for (int i=0; i<nSamples; i++)
		{
			float *gpp = gp_vals+ i*nGenotypes;
			float *glp = gl_vals+ i*nGenotypes;
			float sum = 0;
			for (int j=0; j<nGenotypes; j++)
			{
				if ( bcf_float_is_vector_end(glp[j]) ) {
					NEVER;// we never expect to truncate the vector and finish early
				}
				if ( bcf_float_is_missing(glp[j]) ) {
					bcf_float_set_missing(gpp[j]);
					continue;
				}else{
					gpp[j]=glp[j];
				}
				gpp[j] = pow(10, gpp[j]);
				sum += gpp[j];
			}
			if ( sum<=0 ) continue;
			for (int j=0; j<nGenotypes; j++)
			{
				if ( bcf_float_is_missing(gpp[j]) ) {
					continue;
				}
				if ( bcf_float_is_vector_end(gpp[j]) ) {
					NEVER;
					// break;
				}

				gpp[j] /= sum;
			}
		}
		ASSERT(0==(bcf_update_format_float(out_hdr,out_bcf_rec,"GP",gp_vals,bcf_tags[GP].n)));
		free(gp_vals);
		gp_vals=NULL;
	}


	float* i16_vals=NULL;
	if(1==args->addI16){

		i16_vals= bcf_tag_alloc<float>(I16,0.0);


		// refnref_n_q13_bases[][]
		// [ref|nonref][forward|reverse]

		// 1   #reference Q13 bases on the forward strand
		i16_vals[0] = refnref_n_q13_bases[0][0];

		// 2   #reference Q13 bases on the reverse strand
		i16_vals[1] = refnref_n_q13_bases[0][1];

		// 3   #non-ref Q13 bases on the forward strand
		i16_vals[2] = refnref_n_q13_bases[1][0];
		// 4   #non-ref Q13 bases on the reverse strand
		i16_vals[3] = refnref_n_q13_bases[1][1];

		// 5   sum of reference base qualities
		i16_vals[4] = refnref_sum_qs[0];

		// 6   sum of squares of reference base qualities
		i16_vals[5] = refnref_sum_qs_sq[0];

		// 7   sum of non-ref base qualities
		i16_vals[6] = refnref_sum_qs[1];

		// 8   sum of squares of non-ref base qualities
		i16_vals[7] = refnref_sum_qs_sq[1];


		// 9   sum of ref mapping qualities
		i16_vals[8] = 0;

		// 10  sum of squares of ref mapping qualities
		i16_vals[9] = 0;

		// 11  sum of non-ref mapping qualities
		i16_vals[10] = 0;

		// 12  sum of squares of non-ref mapping qualities
		i16_vals[11] = 0;

		// 13  sum of tail distance for ref bases
		i16_vals[12] = refnref_sum_taildist[0];

		// 14  sum of squares of tail distance for ref bases
		i16_vals[13] = refnref_sum_taildist_sq[0];

		// 15  sum of tail distance for non-ref bases
		i16_vals[14] = refnref_sum_taildist[1];

		// 16  sum of squares of tail distance for non-ref
		i16_vals[15] = refnref_sum_taildist_sq[1];

		
		ASSERT(0==(bcf_update_info_float(out_hdr, out_bcf_rec, "I16", i16_vals,16)));
		free(i16_vals);
		i16_vals=NULL;
	}



	// [BEGIN ADD QS TAG] ------------------------------------------------------- //
	if(1==args->addQS){
		if(NULL==qs_vals){
			NEVER;
		}else{
			ASSERT(0==(bcf_update_info_float(out_hdr, out_bcf_rec, "QS", qs_vals, bcf_tags[QS].n)));

			free(qs_vals);
			qs_vals=NULL;
		}
	}
	// [END ADD QS TAG] --------------------------------------------------------- //


	free(gl_vals);
	gl_vals=NULL;
	free(dp_vals);
	dp_vals=NULL;

	free(gls);
	gls=NULL;

	return(0);
}

int simulate_record_true_values(bcf_hdr_t *out_hdr,bcf1_t *out_bcf_rec,int nSamples,argStruct* args, const int site_i){

	if (0!=args->errate){
		ERROR("Cannot simulate true values when error rate is defined. Please set error rate to 0 and rerun.");
	}

	// without error, only 3 possible genotypes
	// input:
	//	REF	ALT
	//	0	1
	//
	// gts:
	// 00 01 11
	int nGenotypes=3;

	int32_t ngt_arr=0;
	int ngt=bcf_get_genotypes(out_hdr, out_bcf_rec, &gt_arr, &ngt_arr);
	if ( ngt<=0 ){
		ERROR("Could not find GT tag.");
	}
	if(2!=ngt/nSamples){
		ERROR("Ploidy %d is not supported.\n",ngt/nSamples);
	}


	float *gl_vals=NULL;
	gl_vals= bcf_tag_alloc<float>(GL,-0.0,nSamples*nGenotypes);

	int32_t *pl_vals=NULL;
	if(1==args->addPL){
		pl_vals = bcf_tag_alloc<int32_t>(PL,0, nSamples*nGenotypes);
	}

	float* gp_vals=NULL;
	if(1==args->addGP){
		gp_vals=bcf_tag_alloc<float>(GP,0.0, nSamples*nGenotypes);
	}


	for (int sample_i=0; sample_i<nSamples; sample_i++) {

		if(args->printBaseCounts==1){
			NEVER;
		}

		int32_t *gt_ptr = gt_arr + sample_i*SIM_PLOIDY;

		int bin_gts[2] = {0};


		//binary input genotypes from simulated input
		for (int i=0; i<SIM_PLOIDY;i++){
			bin_gts[i]= bcf_gt_allele(gt_ptr[i]);
			// use bit shifting to check if bin_gt[n] is 0 or 1
			if( (bin_gts[i] >> 1) != 0){
				fprintf(stderr,"ERROR:\n\nbin_gts[%d]: Genotype %d not supported.\n",i,bin_gts[i]);
				exit(1);
			}
		}

		// most likely value should be at the genotype combination for true genotype
		int true_gt_idx=bin_gts[0]+bin_gts[1];
		for (int i=0;i<nGenotypes;++i){
			if(i==true_gt_idx){
				gl_vals[sample_i * nGenotypes + i]=MAXGL;

				if(1==args->addPL){
					pl_vals[sample_i * nGenotypes + i]=MAXPL;
				}
				if(1==args->addGP){
					gp_vals[sample_i * nGenotypes + i]=MAXGP;
				}

			}else{
				gl_vals[sample_i * nGenotypes + i]=MINGL;

				if(1==args->addPL){
					pl_vals[sample_i * nGenotypes + i]=MINPL;
				}
				if(1==args->addGP){
					gp_vals[sample_i * nGenotypes + i]=MINGP;
				}
			}
		}



	} //end sample loop


	// set REF and ALT alleles
	ASSERT(0==(bcf_update_alleles_str(out_hdr,out_bcf_rec,"A,C")));



	ASSERT(0==(bcf_update_format_float(out_hdr, out_bcf_rec, "GL", gl_vals,bcf_tags[GL].n)));
	free(gl_vals);
	gl_vals=NULL;


	if(1==args->addGP){
		ASSERT(0==(bcf_update_format_float(out_hdr,out_bcf_rec,"GP",gp_vals,bcf_tags[GP].n)));
		free(gp_vals);
		gp_vals=NULL;
	}

	if(1==args->addPL){
		ASSERT(0==(bcf_update_format_int32(out_hdr,out_bcf_rec,"PL",pl_vals,bcf_tags[PL].n)));
		free(pl_vals);
		pl_vals=NULL;
	}

	if(1==args->addI16){
		NEVER;
	}

	if(1==args->addQS){
		NEVER;
	}

	return(0);
}

int simulate_record(bcf_hdr_t *out_hdr,bcf1_t *out_bcf_rec,int nSamples,double* mps_depths, argStruct* args, const int site_i, FILE* out_baseCounts_ff){


	if(-999==args->mps_depth){
		ASSERT(NULL==mps_depths);
		ASSERT(NULL==out_baseCounts_ff);
		return simulate_record_true_values(out_hdr,out_bcf_rec,nSamples,args, site_i);
	}

	return simulate_record_values(out_hdr,out_bcf_rec,nSamples,mps_depths, args, site_i, out_baseCounts_ff);
}


int main(int argc, char **argv) {


	if(argc==1){
		help_page();
		return 0;
	}
	argStruct *args=args_get(--argc,++argv);

	if(args!=NULL){

		char *in_fn=args->in_fn;
		char *out_fp=args->out_fp;
		char* in_mps_depths=args->in_mps_depths;
		double *mps_depths=NULL;
		int pos0=args->pos0;

		FILE *arg_ff=openFile(out_fp,".arg");

		FILE *out_baseCounts_ff=NULL;
		if(args->printBaseCounts==1){
			out_baseCounts_ff=openFile(out_fp,".baseCounts.tsv");
			fprintf(out_baseCounts_ff,"site\tind\tA\tC\tG\tT\n");
		}

		char *COMMAND=NULL;
		char *SOURCE_TAG;

		if(-999==args->mps_depth){
			ASSERT(asprintf(&COMMAND,"vcfgl --input %s --output %s --output-mode %s --error-rate %f --depth inf --depths-file %s -pos0 %d -seed %d -explode %d -printBaseCounts %d -addGP %d -addPL %d -addI16 %d -addQS %d\n",args->in_fn,args->out_fp,args->output_mode,args->errate,args->in_mps_depths,args->pos0,args->seed,args->explode,args->printBaseCounts,args->addGP,args->addPL,args->addI16,args->addQS)>0);
			fprintf(stderr,"\n%s",COMMAND);
			fprintf(arg_ff,"\n%s",COMMAND);
			ASSERT(asprintf(&SOURCE_TAG, "##source=%s", COMMAND)>0);
		}else{
			ASSERT(asprintf(&COMMAND,"vcfgl --input %s --output %s --output-mode %s --error-rate %f --depth %f --depths-file %s -pos0 %d -seed %d -explode %d -printBaseCounts %d -addGP %d -addPL %d -addI16 %d -addQS %d\n",args->in_fn,args->out_fp,args->output_mode,args->errate,args->mps_depth,args->in_mps_depths,args->pos0,args->seed,args->explode,args->printBaseCounts,args->addGP,args->addPL,args->addI16,args->addQS)>0);
			fprintf(stderr,"\n%s",COMMAND);
			fprintf(arg_ff,"\n%s",COMMAND);
			ASSERT(asprintf(&SOURCE_TAG, "##source=%s", COMMAND)>0);

		}
		free(COMMAND);
		COMMAND=NULL;



		vcfFile * in_ff = bcf_open(in_fn, "r");
		vcfFile * out_ff=NULL;

		char *OUT_EXT=NULL;
		char *out_fn=NULL;
		switch (*args->output_mode){
			case 'v':
				fprintf(stderr,"\nOutput is VCF file\n");
				OUT_EXT=strdup(".vcf");
				out_fn = (char*)malloc(strlen(out_fp)+strlen(OUT_EXT)+1);
				strcpy(out_fn,out_fp);
				strcat(out_fn,OUT_EXT);
				out_ff = bcf_open(out_fn, "w");
				break;
			case 'b':
				fprintf(stderr,"\nOutput is BCF file\n");
				OUT_EXT=strdup(".bcf");
				out_fn = (char*)malloc(strlen(out_fp)+strlen(OUT_EXT)+1);
				strcpy(out_fn,out_fp);
				strcat(out_fn,OUT_EXT);
				out_ff = bcf_open(out_fn, "wb");
				break;
			case 'z':
				fprintf(stderr,"\nOutput is compressed VCF file\n");
				OUT_EXT=strdup(".vcf.gz");
				out_fn = (char*)malloc(strlen(out_fp)+strlen(OUT_EXT)+1);
				strcpy(out_fn,out_fp);
				strcat(out_fn,OUT_EXT);
				out_ff = bcf_open(out_fn, "wz");
				break;
			case 'u':
				fprintf(stderr,"\nOutput is uncompressed BCF file\n");
				OUT_EXT=strdup(".bcf");
				out_fn = (char*)malloc(strlen(out_fp)+strlen(OUT_EXT)+1);
				strcpy(out_fn,out_fp);
				strcat(out_fn,OUT_EXT);
				out_ff = bcf_open(out_fn, "wbu");
				break;
		}
		free(OUT_EXT);

		if (in_ff == NULL) {
			return 1;
		}

		if (bcf == 0) {
			return 1; 
		}


		bcf_hdr_t *hdr = bcf_hdr_read(in_ff);

		bcf_hdr_t *out_hdr = bcf_hdr_dup(hdr);
		bcf_hdr_merge(out_hdr,hdr);


		char *DATE_TAG;
		char *DATETIME=get_time();
		fprintf(stderr,"\n%s\n",DATETIME);

		ASSERT(asprintf(&DATE_TAG, "##fileDate=%s", DATETIME)>0);
		ASSERT(0==bcf_hdr_append(out_hdr, DATE_TAG));
		free(DATE_TAG);
		DATE_TAG=NULL;

		ASSERT(0==bcf_hdr_append(out_hdr, SOURCE_TAG));
		free(SOURCE_TAG);
		SOURCE_TAG=NULL;

		char *SOURCE_VERSION_TAG;
		ASSERT(asprintf(&SOURCE_VERSION_TAG, "##source=vcfgl version: %s",VCFGL_VERSION)>0);
		ASSERT(0==bcf_hdr_append(out_hdr, SOURCE_VERSION_TAG));
		free(SOURCE_VERSION_TAG);
		SOURCE_VERSION_TAG=NULL;


		int nSamples=bcf_hdr_nsamples(hdr);

		ASSERT(0==bcf_hdr_append(out_hdr,bcf_tags[GL].hdr));


		if(1==args->addGP){
			ASSERT(0==bcf_hdr_append(out_hdr,bcf_tags[GP].hdr));
		}

		if(1==args->addPL){
			ASSERT(0==bcf_hdr_append(out_hdr,bcf_tags[PL].hdr));
		}

		if(-999!=args->mps_depth){
			bcf_tag_set_size(DP, nSamples);
			ASSERT(0==bcf_hdr_append(out_hdr,bcf_tags[DP].hdr));

			if(1==args->addI16){
				ASSERT(0==bcf_hdr_append(out_hdr,bcf_tags[I16].hdr));
			}

			if(1==args->addQS){
				ASSERT(0==bcf_hdr_append(out_hdr,bcf_tags[QS].hdr));
			}
		}


		ASSERT(bcf_hdr_write(out_ff,out_hdr)==0);

		bcf1_t *bcf = bcf_init();
		int nSites=0;


		if(pos0){
			fprintf(stderr, "\n -pos0=%d ; This means input VCF's positions are 0 based, and will shift coordinate system with +1\n", pos0);
		}







		if(in_mps_depths!=NULL){
			mps_depths=read_depthsFile(in_mps_depths, nSamples);

			fprintf(stderr, "\n");
			for (int sample_i=0; sample_i<nSamples; sample_i++) {
				fprintf(stderr, "Individual %d mean per-site depth is set to %f\n", sample_i,mps_depths[sample_i]);
			}
		}


		fprintf(stderr, "\nReading file:\t\"%s\"\n", in_fn);
		fprintf(stderr, "Number of samples: %i\n",nSamples); 
		fprintf(stderr,	"Number of contigs: %d\n",hdr->n[BCF_DT_CTG]);


		bcf1_t *out_bcf_rec=bcf_init();
		bcf1_t *blank = bcf_init();

		while (bcf_read(in_ff, hdr, bcf) == 0) {
			//copy next record with data into out_bcf_rec
			out_bcf_rec=bcf_copy(out_bcf_rec,bcf);
			if(args->explode==0){
				ASSERT(0==simulate_record(out_hdr,out_bcf_rec,nSamples,mps_depths,args,nSites,out_baseCounts_ff));
				if(out_bcf_rec->pos==-1){
					if(pos0==0){
						fprintf(stderr,"\n[ERROR]: Input file coordinates start from 0; but -pos0 is not set to 1. Please run again with -pos0 1.\n\n");
						exit(1);
					}
				}
				out_bcf_rec->pos += pos0;
				ASSERT(0==bcf_write(out_ff, out_hdr, out_bcf_rec));
				nSites++;
			}else{
				//ensure that we have a empty blank record that we can modify
				blank = bcf_copy(blank,out_bcf_rec);

				if(out_bcf_rec->pos==-1){
					if(pos0==0){
						fprintf(stderr,"\n[ERROR]: Input file coordinates start from 0; but -pos0 is not set to 1. Please run again with -pos0 1.\n\n");
						exit(1);
					}
				}

				while(1){//this block should run for every site with missing/nodata
						 //  fprintf(stderr,"\t\t-> out_bcf_recpos: %d nSites: %d\n",out_bcf_rec->pos,nSites);
					if(nSites==out_bcf_rec->pos+pos0){
						// fprintf(stderr,"now breaking\n");
						break;
					}

					setblank(blank,out_bcf_rec,hdr);
					ASSERT(0==simulate_record(out_hdr,blank,nSamples,mps_depths,args,nSites,out_baseCounts_ff));
					blank->pos = nSites;
					// fprintf(stderr,"blank->pos: %d\n",blank->pos+pos0);
					if(bcf_write(out_ff, out_hdr, blank)!=0){
						fprintf(stderr,"Error: Failed to write\n");
						exit(1);
					}
					nSites++;

				}
				// fprintf(stderr,"After loop that fills in missing data will print out: %d\n",out_bcf_rec->pos+pos0+1);
				ASSERT(0==simulate_record(out_hdr,out_bcf_rec,nSamples,mps_depths,args,nSites,out_baseCounts_ff));
				out_bcf_rec->pos += pos0;
				if(bcf_write(out_ff, out_hdr, out_bcf_rec)!=0){
					fprintf(stderr,"Error: Failed to write\n");
					exit(1);
				}
				nSites++;
			}
		}

		if(args->explode==1){
			bcf_idpair_t *ctg = hdr->id[BCF_DT_CTG];
			int contigsize = ctg[out_bcf_rec->rid].val->info[0];

			while(nSites<contigsize){
				setblank(blank,out_bcf_rec,hdr);
				ASSERT(0==simulate_record(out_hdr,blank,nSamples,mps_depths,args,nSites,out_baseCounts_ff));
				blank->pos = nSites;
				if(bcf_write(out_ff, out_hdr, blank)!=0){
					fprintf(stderr,"Error: Failed to write\n");
					exit(1);
				}
				nSites++;
			}
		}


		fprintf(stderr, "Total number of sites: %i\n", nSites);

		bcf_hdr_destroy(hdr);
		bcf_destroy(bcf);
		bcf_hdr_destroy(out_hdr);
		bcf_destroy(out_bcf_rec);
		bcf_destroy(blank);

		ASSERT(0==bcf_close(in_ff));
		ASSERT(0==bcf_close(out_ff));
		ASSERT(0==fclose(arg_ff));


		if(args->printBaseCounts==1){
			ASSERT(0==fclose(out_baseCounts_ff));
			fprintf(stderr,"\nDumping baseCounts file to %s.baseCounts.tsv\n",out_fp);
		}

		args_destroy(args);

		free(out_fn);
		out_fn=NULL;

		free(gt_arr);
		gt_arr=NULL;

	}

	return 0;

}

