#include "io.h"

FILE *getFILE(const char *fname, const char *mode)
{
	FILE *fp;
	if (NULL == (fp = fopen(fname, mode)))
	{
		fprintf(stderr, "[%s:%s()]\t->Error opening FILE handle for file:%s exiting\n", __FILE__, __FUNCTION__, fname);
		exit(0);
	}
	return fp;
}

FILE *openFILE(const char *a, const char *b)
{
	char *c = (char *)malloc(strlen(a) + strlen(b) + 1);
	strcpy(c, a);
	strcat(c, b);
	fprintf(stderr, "\t-> Opening output file for writing: %s\n", c);
	FILE *fp = getFILE(c, "w");
	free(c);
	return fp;
}

void help_page()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "vcfgl [version: %s] [build: %s %s] [htslib: %s]\n", VCFGL_VERSION, __DATE__, __TIME__, hts_version());
	fprintf(stderr, "\n");
	fprintf(stderr, "\n");

	fprintf(stderr, "Usage:\tvcfgl -i <input> [options]\n\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Options:\n");

	fprintf(stderr, "\t-v/--verbose\t\tVerbose (default:0=disabled)\n");
	fprintf(stderr, "\t-@/--threads\t\tNumber of threads (default:1)\n");

	fprintf(stderr, "\n");

	fprintf(stderr, "    -->\tInput/Output\n");

	fprintf(stderr, "\t-i/--input\t\tInput file (required)\n");
	fprintf(stderr, "\t-o/--output\t\tOutput file prefix (default:output)\n");
	fprintf(stderr, "\t-O/--output-mode\tOutput mode (default:b)\n");
	fprintf(stderr, "\t\t\t\tv\tVCF file\n");
	fprintf(stderr, "\t\t\t\tb\tBCF file\n");
	fprintf(stderr, "\t\t\t\tz\tCompressed VCF file (vcf.gz)\n");
	fprintf(stderr, "\t\t\t\tu\tUncompressed BCF file\n");
	fprintf(stderr, "\n");

	fprintf(stderr, "    -->\tSimulation parameters\n");
	fprintf(stderr, "\t-d/--depth\t\tMean per-site read depth (default:1.0)\n");
	fprintf(stderr, "\t\t\t\tUse `--depth inf` to set the simulated values to known true variables\n");
	fprintf(stderr, "\t-df/--depths-file\tFile containing mean per-site read depth for each sample (conflicts with -d)\n");
	fprintf(stderr, "\t-e/--error-rate\t\tError rate (default:0.01)\n");
	fprintf(stderr, "\t--error-qs\t\tShould the program sample errors?\n");
	fprintf(stderr, "\t\t\t\t0: no\n");
	fprintf(stderr, "\t\t\t\t1: sample error probabilities from beta distribution\n");
	fprintf(stderr, "\t--beta-variance\t\tVariance of the beta distribution\n");

	fprintf(stderr, "\t--seed\t\t\tRandom seed used to initialize the random number generator\n");
	fprintf(stderr, "\n");

	fprintf(stderr, "\t--pos0\t\t\tAre the input coordinates are 0-based? (default:0=no)");
	fprintf(stderr, "\n\t\t\t\tIf input cordinates are 0 based, use --pos0 1 to shift positions by +1\n");
	fprintf(stderr, "\t--rm-invar-sites\tRemove sites where only one base for observed in total (default:0=disabled)\n");
	// fprintf(stderr, "\t--trim-alt-alleles\tTrim ALT alleles not sampled during base sampling (default:0=disabled)\n");
	fprintf(stderr, "\t--platform\t\tSimulate base qualities for a specific sequencing platform (default:0=disabled)\n");
	fprintf(stderr, "\t\t\t\t1\tNovaSeq 6000. The qualities are binned into four possible quality values: 2, 12, 23 and 37.\n");

	fprintf(stderr, "\n");
	fprintf(stderr, "Commands:\n");

	fprintf(stderr, "    -->\tSpecify which fields to simulate\n");
	fprintf(stderr, "\t-addGL\t\t\tAdd GL field (default:1=enabled)\n");
	fprintf(stderr, "\t-addGP\t\t\tAdd GP field (default:0=disabled)\n");
	fprintf(stderr, "\t-addPL\t\t\tAdd PL field (default:0=disabled)\n");
	fprintf(stderr, "\t-addI16\t\t\tAdd I16 field (default:0=disabled)\n");
	fprintf(stderr, "\t-addQS\t\t\tAdd QS field (default:0=disabled)\n");
	fprintf(stderr, "\t-addFormatDP\t\tAdd FORMAT/DP field (default:1=enabled)\n");
	fprintf(stderr, "\t-addFormatAD\t\tAdd FORMAT/AD field (default:0=disabled)\n");
	fprintf(stderr, "\t-addInfoDP\t\tAdd INFO/DP field (default:0=disabled)\n");
	// fprintf(stderr, "\t-addFormatADF\t\tAdd FORMAT/ADF field (default:0=disabled)\n");
	// fprintf(stderr, "\t-addFormatADR\t\tAdd FORMAT/ADR field (default:0=disabled)\n");

	fprintf(stderr, "\t-explode\t\tAlso simulate sites that were not observed in the input file (default:0=disabled)\n");

	fprintf(stderr, "\n");
}

argStruct *args_init()
{

	args = (argStruct *)calloc(1, sizeof(argStruct));


	args->n_threads = 1;

	args->verbose = 0;

	args->out_fnp = NULL;
	args->in_fn = NULL;

	args->betaSampler = NULL;

	args->datetime = NULL;
	args->command = NULL;

	args->error_qs = 0;
	args->beta_variance = 0.0;

	// for setting the quality score with a fixed given error rate
	args->error_rate_q = -1;

	args->pos0 = 0;
	args->trimAlts = 0;
	args->useUnknownAllele =0;
	args->platform = 0;

	args->rmInvarSites = 0;

	args->addGL = 1;
	args->addGP = 0;
	args->addPL = 0;
	args->addI16 = 0;
	args->addQS = 0;
	args->addFormatDP = 1;
	args->addInfoDP= 0;
	args->addFormatAD = 0;
	args->addFormatADF = 0; //TODO
	args->addFormatADR = 0; //TODO
	args->addInfoAD = 0; 
	args->addInfoADF = 0; //TODO
	args->addInfoADR = 0; //TODO

	args->mps_depth = 1.0;
	args->mps_depths_fn = NULL;
	args->error_rate = 0.01;
	args->seed = -1;
	args->explode = 0;
	return args;
}

void args_destroy(argStruct *args)
{

	free(args->in_fn);
	args->in_fn = NULL;

	free(args->out_fnp);
	args->out_fnp = NULL;

	free(args->out_fn);
	args->out_fn = NULL;

	free(args->output_mode);
	args->output_mode = NULL;

	free(args->output_mode_str);
	args->output_mode_str = NULL;

	free(args->datetime);
	args->datetime = NULL;

	free(args->command);
	args->command = NULL;

	if (NULL != args->mps_depths_fn)
	{
		free(args->mps_depths_fn);
		args->mps_depths_fn = NULL;
	}

	delete args->betaSampler;

	free(args);
	args = NULL;
}

char *get_time()
{
	time_t current_time;
	struct tm *local_time;
	current_time = time(NULL);
	local_time = localtime(&current_time);
	return (asctime(local_time));
}

argStruct *args_get(int argc, char **argv)
{

	args = args_init();

	if (0 == argc)
	{
		help_page();
		exit(0);
	}

	while (*argv)
	{

		char *arv = *argv;
		char *val = *(++argv);
		if ((strcmp("-h", arv) == 0) || (strcmp("--help", arv) == 0))
		{
			help_page();
			exit(0);
		}

		else if ((strcmp("--verbose", arv) == 0) || (strcmp("-v", arv) == 0))
		{
			args->verbose = atoi(val);
		}

		else if ((strcmp("--threads", arv) == 0) || (strcmp("-@", arv) == 0))
		{
			args->n_threads = atoi(val);
		}

		else if ((strcmp("-i", arv) == 0) || (strcmp("--input", arv) == 0))
		{
			args->in_fn = strdup(val);
		}

		else if ((strcmp("-o", arv) == 0) || (strcmp("--output", arv) == 0))
		{
			args->out_fnp = strdup(val);
		}

		else if ((strcmp("-O", arv) == 0) || (strcmp("--output-mode", arv) == 0))
		{
			args->output_mode = strdup(val);
		}

		else if ((strcmp("-e", arv) == 0) || (strcmp("--error-rate", arv) == 0))
		{
			args->error_rate = atof(val);
		}

		else if ((strcmp("-d", arv) == 0) || (strcmp("--depth", arv) == 0))
		{
			char *tmp = strdup(val);
			bool is_number = true;
			int i = 0;
			if (tmp[0] == '-')
				i = 1;
			for (; tmp[i] != 0; ++i)
			{
				if (!isdigit(tmp[i]))
				{
					if (tmp[i] == '.')
						continue;
					is_number = false;
					if (strcasecmp("inf", tmp + i) == 0)
					{
						args->mps_depth = -999;
					}
				}
			}
			free(tmp);
			tmp = NULL;
			if (is_number)
			{
				args->mps_depth = atof(val);
			}
			else if (-999 != args->mps_depth)
			{
				ERROR("--depth is set to unknown value %s", val);
			}
		}

		else if ((strcmp("-df", arv) == 0) || (strcmp("--depths-file", arv) == 0))
		{
			args->mps_depths_fn = strdup(val);
			args->mps_depth = -1;
		}

		else if (strcasecmp("--error-qs", arv) == 0)
		{
			args->error_qs = atoi(val);
		}
		else if (strcasecmp("--beta-variance", arv) == 0)
		{
			args->beta_variance = atof(val);
		}

		else if (strcasecmp("--pos0", arv) == 0)
			args->pos0 = atoi(val);

		else if ((strcmp("-s", arv) == 0) || (strcmp("--seed", arv) == 0))
		{
			args->seed = atoi(val);
		}

		else if (strcmp("-explode", arv) == 0)
			args->explode = atoi(val);
		else if (strcasecmp("-addGL", arv) == 0)
			args->addGL = atoi(val);
		else if (strcasecmp("-addGP", arv) == 0)
			args->addGP = atoi(val);
		else if (strcasecmp("-addPL", arv) == 0)
			args->addPL = atoi(val);
		else if (strcasecmp("-addI16", arv) == 0)
			args->addI16 = atoi(val);
		else if (strcasecmp("-addQS", arv) == 0)
			args->addQS = atoi(val);
		else if (strcasecmp("-addFormatDP", arv) == 0)
			args->addFormatDP = atoi(val);
		else if (strcasecmp("-addFormatAD", arv) == 0)
			args->addFormatAD = atoi(val);
		else if (strcasecmp("-addFormatADF", arv) == 0)
			args->addFormatADF = atoi(val);
		else if (strcasecmp("-addFormatADR", arv) == 0)
			args->addFormatADR = atoi(val);
		else if (strcasecmp("-addInfoDP", arv) == 0)
			args->addInfoDP = atoi(val);
		else if (strcasecmp("-addInfoAD", arv) == 0)
			args->addInfoAD = atoi(val);
		else if (strcasecmp("-addInfoADF", arv) == 0)
			args->addInfoADF = atoi(val);
		else if (strcasecmp("-addInfoADR", arv) == 0)
			args->addInfoADR = atoi(val);

		else if (strcasecmp("--rm-invar-sites", arv) == 0)
			args->rmInvarSites = atoi(val);

		else if (strcasecmp("--trim-alt-alleles", arv) == 0)
			args->trimAlts = atoi(val);

		else if (strcasecmp("--use-unknown-allele", arv) == 0)
			args->useUnknownAllele = atoi(val);
		

		else if (strcasecmp("--platform", arv) == 0)
		{
			args->platform = atoi(val);
		}

		else
		{
			fprintf(stderr, "Unknown arg:%s\n", arv);
			free(args);
			return 0;
		}
		++argv;
	}


	
	if(args->rmInvarSites){

	//TODO
	//
		NEVER;
	}

	if (1 == args->addI16)
	{
		WARN("I16 tag simulation is still under development. Please use with caution.");
		WARN("-addI16 1 is used. Will set mapping quality-related I16 values to 0.");
	}

	if (1==args->trimAlts)
	{
		NEVER;
	}

	if(1==args->explode){
		// NEVER;
	}

	if(1!=args->n_threads){
		if(0==args->n_threads){
			VWARN("--threads 0 implies --threads 1. Setting number of threads to 1 (no multithreading).");
			args->n_threads=1;
		}
	}

	if (args->seed == -1)
	{
		int rseed = time(NULL);
		fprintf(stderr, "\n-> No seed was given. Setting the random seed to the randomly chosen value: %d\n", rseed);
		args->seed = rseed;
	}
	srand48(args->seed);

	if (args->in_fn == NULL)
	{
		ERROR("Input file is not specified. Please use -i/--input option to specify the input file.");
	}

	if (NULL == args->output_mode)
	{
		args->output_mode = strdup("b");
	}

	if (NULL == args->out_fnp)
	{
		args->out_fnp = strdup("output");
	}

	args->datetime = strdup(get_time());

	// -depth inf
	if (-999 == args->mps_depth)
	{
		
		
		if(args->addFormatDP){
			ERROR("(-addFormatDP 1) FORMAT/DP tag cannot be added when --depth inf is set.");
		}
		
		if(args->addInfoDP){
			ERROR("(-addInfoDP 1) INFO/DP tag cannot be added when --depth inf is set.");
		}


		if (1 == args->addI16)
		{
			ERROR("(-addI16 1) I16 tag cannot be added when --depth inf is set.");

		}

		if (1 == args->addQS)
		{
			ERROR("(-addQS 1) QS tag cannot be added when --depth inf is set.");
			//TODO maybe add a very very high qs val instead?
		}

		if (args->addFormatAD){
			ERROR("(-addFormatAD 1) FORMAT/AD tag cannot be added when --depth inf is set.");
		}
		if (args->addFormatADF){
			ERROR("(-addFormatADF 1) FORMAT/ADF tag cannot be added when --depth inf is set.");
		}
		if (args->addFormatADR){
			ERROR("(-addFormatADR 1) FORMAT/ADR tag cannot be added when --depth inf is set.");
		}

		if (args->addInfoAD){
			ERROR("(-addInfoAD 1) INFO/AD tag cannot be added when --depth inf is set.");
		}
		if (args->addInfoADF){
			ERROR("(-addInfoADF 1) INFO/ADF tag cannot be added when --depth inf is set.");
		}
		if (args->addInfoADR){
			ERROR("(-addInfoADR 1) INFO/ADR tag cannot be added when --depth inf is set.");
		}
		

		if (0 != args->error_rate)
		{
			ERROR("Cannot simulate true values (--depth inf) with error rate (--error-rate) set to %f. Please set --error-rate to 0 and rerun.", args->error_rate);
		}

	}

	if (1 == args->pos0)
	{
		fprintf(stderr, "\n --pos0=%d ; This means input VCF's positions are 0 based, and will shift coordinate system with +1\n", args->pos0);
	}

	if (1 == args->error_qs)
	{

		if (0 == args->error_rate)
		{
			ERROR("--error-qs %d requires --error-rate to be set. Please set --error-rate and rerun.", args->error_qs);
		}

		if (args->beta_variance < 0)
		{
			ERROR("--error-qs %d requires --beta-variance to be set to a positive value. Please set --beta-variance and rerun.", args->error_qs);
		}else if(0==args->beta_variance){
			WARN("--error-qs %d, but --beta-variance is set to 0. This is equivalent to --error-qs 0. Setting --error-qs to 0.", args->error_qs);
			args->error_qs=0;
		}else{

			args->betaSampler = new BetaSampler(args->error_rate, args->beta_variance, args->seed);
		}

	}


	if ((0 != args->beta_variance)&&(0==args->error_qs)){
			ERROR("--beta-variance %d requires --error-qs 1.", args->error_qs);
	}

	// TODO but it should be possible to set the error rate to 0 without --depth inf
	// error rate can only be 0 with --depth inf
	// if (-999 != args->mps_depth && 0 == args->error_rate)
	// {
	// ERROR("--error-rate %f requires --depth inf. Please set --depth inf and rerun.", args->error_rate);
	// }

	char depth_val[100];
	if (-999 == args->mps_depth)
	{
		sprintf(depth_val, "%s", "inf");
	}
	else if (-1 == args->mps_depth)
	{
		// depths file is defined
		sprintf(depth_val, "%s", "depths_file");
	}
	else
	{
		sprintf(depth_val, "%f", args->mps_depth);
	}

	ASSERT(asprintf(&args->command, "vcfgl --input %s --output %s --output-mode %s --error-rate %f --depth %s --depths-file %s --pos0 %d --seed %d --error-qs %d --beta-variance %e --rm-invar-sites %d --trim-alt-alleles %d --platform %d -explode %d -addGL %d -addGP %d -addPL %d -addI16 %d -addQS %d  -addFormatDP %d -addFormatAD %d -addFormatADF %d -addFormatADR %d -addInfoDP %d -addInfoAD %d -addInfoADF %d -addInfoADR %d", args->in_fn, args->out_fnp, args->output_mode, args->error_rate, depth_val, args->mps_depths_fn, args->pos0, args->seed, args->error_qs, args->beta_variance, args->rmInvarSites, args->trimAlts, args->platform, args->explode, args->addGL, args->addGP, args->addPL, args->addI16, args->addQS, args->addFormatDP, args->addFormatAD, args->addFormatADF, args->addFormatADR, args->addInfoDP, args->addInfoAD, args->addInfoADF, args->addInfoADR) > 0);


	args->out_fn = NULL;
	args->output_mode_str=NULL;


	switch (*args->output_mode)
	{
	case 'v':
		fprintf(stderr, "\nOutput is VCF file (.vcf)\n");
		if(args->n_threads>1){
			ERROR("Multithreading is not supported for VCF output. Please set --threads 1 and rerun.");
		}
		args->out_fn = (char *)malloc(strlen(args->out_fnp) + strlen(".vcf") + 1);
		strcpy(args->out_fn, args->out_fnp);
		strcat(args->out_fn, ".vcf");
		args->output_mode_str=strdup("w");
		break;
	case 'b':
		fprintf(stderr, "\nOutput is BCF file (.bcf)\n");
		args->out_fn = (char *)malloc(strlen(args->out_fnp) + strlen(".bcf") + 1);
		strcpy(args->out_fn, args->out_fnp);
		strcat(args->out_fn, ".bcf");
		args->output_mode_str=strdup("wb");
		break;
	case 'z':
		if(args->n_threads>1){
			ERROR("Multithreading is not supported for VCF output. Please set --threads 1 and rerun.");
		}
		fprintf(stderr, "\nOutput is compressed VCF file (.vcf.gz)\n");
		args->out_fn = (char *)malloc(strlen(args->out_fnp) + strlen(".vcf.gz") + 1);
		strcpy(args->out_fn, args->out_fnp);
		strcat(args->out_fn, ".vcf.gz");
		args->output_mode_str=strdup("wz");
		break;
	case 'u':
		fprintf(stderr, "\nOutput is uncompressed BCF file (.bcf)\n");
		args->out_fn = (char *)malloc(strlen(args->out_fnp) + strlen(".bcf") + 1);
		strcpy(args->out_fn, args->out_fnp);
		strcat(args->out_fn, ".bcf");
		args->output_mode_str=strdup("wbu");
		break;
	}

	return (args);
}

size_t fsize(const char *fname)
{
	struct stat st;
	stat(fname, &st);
	return st.st_size;
}

// modified from msToGlf.c
double *read_depthsFile(const char *fname, int len)
{
	if (args->verbose > 0)
	{
		fprintf(stderr, "Reading depths file: %s for %d samples\n", fname, len);
	}

	FILE *fp = NULL;
	if ((fp = fopen(fname, "r")) == NULL)
	{
		fprintf(stderr, "Problem opening file: %s\n", fname);
		exit(0);
	}

	char buf[1024];
	double *ret = (double *)malloc(len * sizeof(double));
	if (fsize(fname) != fread(buf, sizeof(char), fsize(fname), fp))
	{
		fprintf(stderr, "Problem reading file=%s\n", fname);
		exit(0);
	}
	int posi = 0;

	ret[posi++] = atof(strtok(buf, "\n\t "));
	if (len > 1)
	{
		char *tok;
		while ((tok = strtok(NULL, "\n\t ")))
			ret[posi++] = atof(tok);
	}

	fclose(fp);

	return ret;
}
