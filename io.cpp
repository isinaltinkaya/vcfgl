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

	fprintf(stderr, "Usage: ./vcfgl -i <input> [options]\n\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Options:\n");

	fprintf(stderr, "\t-i/--input\t\tInput file (required)\n");
	fprintf(stderr, "\t-o/--output\t\tOutput file prefix (default:output)\n");
	fprintf(stderr, "\t-O/--output-mode\t\tOutput mode (default:b)\n");
	fprintf(stderr, "\t\t\t\tv\tVCF file\n");
	fprintf(stderr, "\t\t\t\tb\tBCF file\n");
	fprintf(stderr, "\t\t\t\tz\tCompressed VCF file (vcf.gz)\n");
	fprintf(stderr, "\t\t\t\tb\tUncompressed BCF file\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "\t-e/--error-rate\t\tError rate (default:0.01)\n");
	fprintf(stderr, "\t-d/--depth\t\tMean per-site read depth (default:1.0)\n");
	fprintf(stderr, "\t\t\t\tUse `--depth inf` to set the simulated values to known true variables\n");
	fprintf(stderr, "\t-df/--depths-file\tFile containing mean per-site read depth for each sample (conflicts with -d)\n");
	fprintf(stderr, "\t--seed\t\t\tRandom seed used to initialize the random number generator\n");

	fprintf(stderr, "\t--pos0\t\t\tAre the input coordinates are 0-based? (default:0)");
	fprintf(stderr, "\n\t\t\t\tIf input cordinates are 0 based, use --pos0 1 to shift positions by +1\n");
	fprintf(stderr, "\t--trim-alt-alleles\tTrim ALT alleles not observed in simulated bases (default:0=disabled)\n");

	fprintf(stderr, "\n");
	fprintf(stderr, "\t-explode\t\tExplode to unobserved sites in the input file (default:0=disabled)\n");
	fprintf(stderr, "\t-printBaseCounts\tPrint base counts (default:0=disabled)\n");
	fprintf(stderr, "\t-addGP\t\t\tAdd GP field (default:0=disabled)\n");
	fprintf(stderr, "\t-addPL\t\t\tAdd PL field (default:0=disabled)\n");
	fprintf(stderr, "\t-addI16\t\t\tAdd I16 field (default:0=disabled)\n");
	fprintf(stderr, "\t-addQS\t\t\tAdd QS field (default:0=disabled)\n");

	fprintf(stderr, "\n");
}

argStruct *args_init()
{

	args = (argStruct *)calloc(1, sizeof(argStruct));

	args->out_fp = NULL;
	args->in_fn = NULL;

	args->datetime = NULL;
	args->command = NULL;

	args->pos0 = 0;
	args->trimAlts = 0;

	args->addGP = 0;
	args->addPL = 0;
	args->addI16 = 0;
	args->addQS = 0;

	args->mps_depth = 1.0;
	args->in_mps_depths = NULL;
	args->errate = 0.01;
	args->seed = -1;
	args->explode = 0;
	args->printBaseCounts = 0;
	return args;
}

void args_destroy(argStruct *args)
{

	free(args->in_fn);
	args->in_fn = NULL;

	free(args->out_fp);
	args->out_fp = NULL;

	free(args->output_mode);
	args->output_mode = NULL;

	free(args->datetime);
	args->datetime = NULL;

	free(args->command);
	args->command = NULL;

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

		else if ((strcmp("-i", arv) == 0) || (strcmp("--input", arv) == 0))
		{
			args->in_fn = strdup(val);
		}

		else if ((strcmp("-o", arv) == 0) || (strcmp("--output", arv) == 0))
		{
			args->out_fp = strdup(val);
		}

		else if ((strcmp("-O", arv) == 0) || (strcmp("--output-mode", arv) == 0))
		{
			args->output_mode = strdup(val);
		}

		else if ((strcmp("-e", arv) == 0) || (strcmp("--error-rate", arv) == 0))
		{
			args->errate = atof(val);
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
			args->in_mps_depths = strdup(val);
		}

		else if (strcasecmp("--pos0", arv) == 0)
			args->pos0 = atoi(val);

		else if ((strcmp("-s", arv) == 0) || (strcmp("--seed", arv) == 0))
		{
			args->seed = atoi(val);
		}

		else if (strcmp("-explode", arv) == 0)
			args->explode = atoi(val);
		else if (strcasecmp("-printBaseCounts", arv) == 0)
			args->printBaseCounts = atoi(val);
		else if (strcasecmp("-addGP", arv) == 0)
			args->addGP = atoi(val);
		else if (strcasecmp("-addPL", arv) == 0)
			args->addPL = atoi(val);
		else if (strcasecmp("-addI16", arv) == 0)
			args->addI16 = atoi(val);
		else if (strcasecmp("-addQS", arv) == 0)
			args->addQS = atoi(val);

		else if (strcasecmp("--trim-alt-alleles", arv) == 0)
			args->trimAlts = atoi(val);

		else
		{
			fprintf(stderr, "Unknown arg:%s\n", arv);
			free(args);
			return 0;
		}
		++argv;
	}

	if (1 == args->trimAlts && 1 != args->addQS)
	{
		ERROR("--trim-alt-alleles requires -addQS 1.");
	}

	if (1 == args->addI16)
	{
		WARNING("-addI16 1 is used. Will set mapping quality-related I16 values to 0.");
	}

	if (args->seed == -1)
	{
		int rseed = time(NULL);
		fprintf(stderr, "\n-> No seed was given. Setting the random seed to the randomly chosen value: %d", rseed);
		args->seed = rseed;
	}
	srand48(args->seed);

	if (args->in_fn == NULL)
	{
		fprintf(stderr, "Must supply -in\n");
		free(args);
		return 0;
	}

	if (NULL == args->output_mode)
	{
		args->output_mode = strdup("b");
	}

	if (NULL == args->out_fp)
	{
		args->out_fp = strdup("output");
	}

	args->datetime = strdup(get_time());

	if (-999 == args->mps_depth)
	{
		if (0 != args->errate)
		{
			ERROR("Cannot simulate true values when error rate is defined. Please set error rate to 0 and rerun.");
		}
		ASSERT(asprintf(&args->command, "vcfgl --input %s --output %s --output-mode %s --error-rate %f --depth inf --depths-file %s --pos0 %d --seed %d -explode %d -printBaseCounts %d -addGP %d -addPL %d -addI16 %d -addQS %d\n", args->in_fn, args->out_fp, args->output_mode, args->errate, args->in_mps_depths, args->pos0, args->seed, args->explode, args->printBaseCounts, args->addGP, args->addPL, args->addI16, args->addQS) > 0);
	}
	else
	{
		ASSERT(asprintf(&args->command, "vcfgl --input %s --output %s --output-mode %s --error-rate %f --depth %f --depths-file %s --pos0 %d --seed %d -explode %d -printBaseCounts %d -addGP %d -addPL %d -addI16 %d -addQS %d\n", args->in_fn, args->out_fp, args->output_mode, args->errate, args->mps_depth, args->in_mps_depths, args->pos0, args->seed, args->explode, args->printBaseCounts, args->addGP, args->addPL, args->addI16, args->addQS) > 0);
	}

	if (1 == args->pos0)
	{
		fprintf(stderr, "\n --pos0=%d ; This means input VCF's positions are 0 based, and will shift coordinate system with +1\n", args->pos0);
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
	fprintf(stderr, "Reading depths file: %s for %d samples\n", fname, len);

	FILE *fp = NULL;
	if ((fp = fopen(fname, "r")) == NULL)
	{
		fprintf(stderr, "Problem opening file: %s\n", fname);
		exit(0);
	}

	char *buf = (char *)malloc(sizeof(char) * fsize(fname));
	double *ret = (double *)malloc(len * sizeof(double));
	if (fsize(fname) != fread(buf, sizeof(char), fsize(fname), fp))
	{
		fprintf(stderr, "Problem reading file=%s\n", fname);
		exit(0);
	}
	int posi = 0;

	ret[posi++] = atof(strtok(buf, "\n\t "));
	char *tok;
	while ((tok = strtok(NULL, "\n\t ")))
		ret[posi++] = atof(tok);

	return ret;
}

kstring_t *kbuf_init()
{
	kstring_t *kbuf = new kstring_t;
	kbuf->l = 0;
	kbuf->m = 0;
	kbuf->s = NULL;
	return kbuf;
}

void kbuf_destroy(kstring_t *kbuf)
{
	if (NULL != kbuf)
	{
		free(kbuf->s);
		kbuf->s = NULL;
		delete kbuf;
	}
}
