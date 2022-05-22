
#ifndef __ARGUMENTS__
#define __ARGUMENTS_




/*
 * @typedef
 * @abstract params - parameter structure
 *
 * @field *in_fn	pointer to input file name
 * @field *out_fp	pointer to output file prefix [angsdput]
 * @field errate	error rate [0.01]
 * @field isSim		is input simulated data? [0]
 *		if input is simulated VCF from tskit.write_vcf
 *		positions are 0-based; shift coordinate system+1;
 * @field mps_depth	mean per site depth [4]
 * @field seed
 *
 *
 */



typedef struct{

	char **argv;


	char* in_fn;
	char* out_fp;

	double errate;
	double mps_depth;

	int isSim;
	int seed;
	
}argStruct;


argStruct *args_init();

argStruct *args_get(int argc, char **argv);




#endif
