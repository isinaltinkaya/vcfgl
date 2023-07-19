
#ifndef __ARGUMENTS__
#define __ARGUMENTS_




/*
 * Macro:[ASSERT]
 * shortcut to evaluate an expression, works the same way as the C-macro assert
 */
#define ASSERT(expr) if (!(expr)) {fprintf(stderr,"\n\n*******\n[ERROR](%s:%d) %s\n*******\n",__FILE__,__LINE__,#expr);exit(1);}



/*
 * @typedef
 * @abstract params - parameter structure
 *
 * @field *in_fn			pointer to input file name
 * @field *out_fp			pointer to output file prefix [angsdput]
 * @field errate			error rate [0.01]
 * @field mps_depth			mean per site depth [1]
 * @field in_mps_depths		assign depths to individuals from per individual mean per site depth file, one line per individual
 * @field pos0				are positions 0-based? [0]
 *							if 1 is set; input VCF positions are 0-based;
 *							shift coordinate system+1;
 * @field seed
 * @field output_mode		char defining the output file format
 *					
 *							Output modes
 *						b	compressed BCF [default]
 *						u	uncompressed BCF
 *						v	uncompressed  VCF
 *						z	compressed VCF (bgzf compressed)
 *
 * @field explode
 * @field printBaseCounts	should the program print base counts
 */
typedef struct{

	char **argv;


	char* in_fn;
	char* out_fp;

	double errate;
	double mps_depth;
	char* in_mps_depths;

	int pos0;
	int seed;



	char* output_mode;

	int explode;
	int printBaseCounts;
	
}argStruct;


argStruct *args_init();

argStruct *args_get(int argc, char **argv);


double *read_depthsFile(const char* fname,int len);


#endif
