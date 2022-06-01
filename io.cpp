#include "io.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>


#include <stdio.h>

argStruct *args_init(){


	argStruct *args=(argStruct*)calloc(1,sizeof(argStruct));

	args->out_fp=strdup("output");


	args->in_fn=NULL;
	args->pos0=0;

	args->output_mode=strdup("b");

	args->mps_depth=4;
	args->errate=0.01;
	args->seed=-1;
	args->in_fa = NULL;
	args->explode = 0;
	return args;


}




argStruct *args_get(int argc, char **argv){

	argStruct *args = args_init(); 

	while(*argv){


		char *arv=*argv;
		char *val=*(++argv);

		if(strcasecmp("-in",arv)==0) args->in_fn=strdup(val);
		else if(strcasecmp("-out",arv)==0) args->out_fp=strdup(val); 
		else if(strcasecmp("-err",arv)==0) args->errate=atof(val); 
		else if(strcasecmp("-depth",arv)==0) args->mps_depth=atof(val); 
		else if(strcasecmp("-pos0",arv)==0) args->pos0=atoi(val);
		else if(strcasecmp("-seed",arv)==0) args->seed=atoi(val);
		else if(strcasecmp("-mode",arv)==0) args->output_mode=strdup(val);
		else if(strcasecmp("-O",arv)==0) args->output_mode=strdup(val);
		else if(strcasecmp("-in_fa",arv)==0) args->in_fa=strdup(val);
		else if(strcasecmp("-explode",arv)==0) args->explode=atoi(val);

		else{
			fprintf(stderr,"Unknown arg:%s\n",arv);
			free(args);
			return 0;
		}
		++argv; 
	} 


	if (args->seed == -1){
		srand48(time(NULL));
	}else{
		srand48(args->seed);
	}

	if(args->in_fn==NULL){
		fprintf(stderr,"Must supply -in\n");
		free(args);
		return 0;
	}



	return args;

}




