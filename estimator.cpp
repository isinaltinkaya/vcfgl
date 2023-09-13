#include "estimator.h"

void gl_ln(int base, double error_rate, double *like)
{

#if 0
	for(int i=0;i<10;i++)
		fprintf(stderr," %f ",like[i]);
#endif

	// fprintf(stderr,"base=%d\n",base);
	/*0=AA, 1=AC, 2=AG, 3=AT, 4=CC, 5=CG, 6=CT, 7=GG, 8= GT, 9= TT*/
	static double homTrue, het, homFalse;
	static int preCalc = 0;
	if (preCalc == 0)
	{
		homTrue = log(1.0 - error_rate);
		het = log((1.0 - error_rate) / 2.0 + error_rate / 6.0);
		homFalse = log(error_rate / 3.0);
		preCalc = 1;
	}
	// fprintf(stderr,"offs=%d homeTrue=%f\n",lut_acgt_offsets[base][0],homTrue);
	/*0=AA*/
	like[lut_acgt_offsets[base][0]] += homTrue; // homozygotic hit

	/*1=AC, 2=AG, 3=AT*/
	for (int o = 1; o < 4; o++) // heterozygotic hit{
		like[lut_acgt_offsets[base][o]] += het;

	/*4=CC, 5=CG, 6=CT, 7=GG, 8= GT, 9= TT*/
	for (int o = 4; o < 10; o++) // non hit
		like[lut_acgt_offsets[base][o]] += homFalse;
}

void gl_log10(int base, double error_rate, double *like)
{

#if 0
	for(int i=0;i<10;i++)
		fprintf(stderr," %f ",like[i]);
#endif

	// fprintf(stderr,"base=%d\n",base);
	/*0=AA, 1=AC, 2=AG, 3=AT, 4=CC, 5=CG, 6=CT, 7=GG, 8= GT, 9= TT*/
	static double homTrue, het, homFalse;
	static int preCalc = 0;
	if (preCalc == 0)
	{
		homTrue = log10(1.0 - error_rate);
		het = log10((1.0 - error_rate) / 2.0 + error_rate / 6.0);
		homFalse = log10(error_rate / 3.0);
		preCalc = 1;
	}
	// fprintf(stderr,"offs=%d homeTrue=%f\n",lut_acgt_offsets[base][0],homTrue);
	/*0=AA*/
	like[lut_acgt_offsets[base][0]] += homTrue; // homozygotic hit

	/*1=AC, 2=AG, 3=AT*/
	for (int o = 1; o < 4; o++) // heterozygotic hit{
		like[lut_acgt_offsets[base][o]] += het;

	/*4=CC, 5=CG, 6=CT, 7=GG, 8= GT, 9= TT*/
	for (int o = 4; o < 10; o++) // non hit
		like[lut_acgt_offsets[base][o]] += homFalse;
}

void rescale_likelihood_ratio(double *like)
{

	// rescale to likeratios
	double mx = like[0];

	for (int i = 1; i < 10; i++)
	{
		if (like[i] > mx)
		{
			mx = like[i];
		}
	}

	for (int i = 0; i < 10; i++)
	{
		like[i] -= mx;
	}
}

void rescale_likelihood_ratio(float *like)
{

	// rescale to likeratios
	float mx = like[0];

	for (int i = 1; i < 10; i++)
	{
		if (like[i] > mx)
		{
			mx = like[i];
		}
	}

	for (int i = 0; i < 10; i++)
	{
		like[i] -= mx;
	}
}

void rescale_likelihood_ratio(float *like, const int size)
{

	ASSERT(size > 0 && size <= MAX_NGTS);
	float mx = like[0];

	for (int i = 1; i < size; i++)
	{
		if (like[i] > mx)
		{
			mx = like[i];
		}
	}

	for (int i = 0; i < size; i++)
	{
		like[i] -= mx;
	}
}
