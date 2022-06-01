
#ifndef __ARGUMENTS__
#define __ARGUMENTS_




/*
 * @typedef
 * @abstract params - parameter structure
 *
 * @field *in_fn	pointer to input file name
 * @field *out_fp	pointer to output file prefix [angsdput]
 * @field errate	error rate [0.01]
 * @field mps_depth	mean per site depth [4]
 * @field pos0		is positions 0-based? [0]
 *		if 1 is set; input VCF positions are 0-based;
 *		shift coordinate system+1;
 * @field seed
 * @field output_mode	char defining the output file format
 *					
 *					Output modes
 *					b	compressed BCF [default]
 *					u	uncompressed BCF
 *					v	uncompressed  VCF
 *					z	compressed VCF (bgzf compressed)
 *
 */



typedef struct{

	char **argv;


	char* in_fn;
	char* out_fp;

	double errate;
	double mps_depth;

	int pos0;
	int seed;



	char* output_mode;

	char *in_fa;
	int explode;
	
}argStruct;


argStruct *args_init();

argStruct *args_get(int argc, char **argv);




#endif
