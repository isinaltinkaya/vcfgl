#include "estimator.h"

#include <stdio.h>
#include <math.h>

const int offsets[4][10]={
	{0,1,2,3,  4,5,6,7,8,9},//AA,AC,AG,AT,therest
	{4,1,5,6,  0,2,3,7,8,9},//CC,AC,CG,CT,therest
	{7,2,5,8,  0,1,3,4,6,9},//GG,AG,CG,GT,therest
	{9,3,6,8,  0,1,2,4,5,7},//TT,AT,CT,GT,therest
};



void gl_log10(int base, double errate, double *like){

#if 0
	for(int i=0;i<10;i++)
		fprintf(stderr," %f ",like[i]);
#endif

	// fprintf(stderr,"base=%d\n",base);
	/*0=AA, 1=AC, 2=AG, 3=AT, 4=CC, 5=CG, 6=CT, 7=GG, 8= GT, 9= TT*/
	static double homTrue,het,homFalse;
	static int preCalc=0;
	if(preCalc==0){
		homTrue = log10(1.0-errate);
		het = log10((1.0-errate)/2.0+errate/6.0);
		homFalse = log10(errate/3.0);
		preCalc = 1;
	}
	//fprintf(stderr,"offs=%d homeTrue=%f\n",offsets[base][0],homTrue);
	/*0=AA*/
	like[offsets[base][0]] += homTrue; //homozygotic hit

	/*1=AC, 2=AG, 3=AT*/
	for(int o=1;o<4;o++)//heterozygotic hit{
		like[offsets[base][o]] += het;

	/*4=CC, 5=CG, 6=CT, 7=GG, 8= GT, 9= TT*/
	for(int o=4;o<10;o++)//non hit
		like[offsets[base][o]] += homFalse;

}


void rescale_likelihood_ratio(double *like){

	//rescale to likeratios
	double mx = like[0];

	for(int i=1;i<10;i++){
		if(like[i]>mx){
			mx=like[i];
		}
	}

	for(int i=0;i<10;i++){
		like[i] -= mx;
	}

}
