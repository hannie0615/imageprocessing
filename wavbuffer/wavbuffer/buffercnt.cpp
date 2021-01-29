#include "spenhbuff.h"
#include "enh_const.h"
#include <string.h>
#include <stdio.h>

// initialization for buffer processing
int initbuff(WAVBUFF *spdata, CAOV *cvdata)
{
	int i, k;	

	// initialization of spdata
	spdata->num_sampling_rate = DEFAULT_SAMPLING_RATE;
	spdata->num_bits_per_sample = 16;
	spdata->bufsize = BUFFSIZE;

	for (i=0; i<2*BUFFSIZE; i++)
	//for (i=0; i<BUFFSIZE; i=i+2)
	{		
		spdata->bufdata[i] = 0;
		spdata->prev_bufdata[i] = 0;	
		spdata->enhdata[i] = 0;		
	}

	for (i=0; i<BUFFSIZE; i++)	
	{	
		spdata->nextenh[i] = 0;		
	}

	// initialization of cvdata
	for(k=0; k < FRSTEP; k++)
	{
		cvdata->xold[k] = 0.0;
	}

	for(k=0; k < FFTSIZE; k++)
	{
		cvdata->noise_ps[k] = 0.0;
		cvdata->prev_ph[k] = 0.5;
		cvdata->prev_specenh[k] = 0.0;
		cvdata->pwsp[k] = 0.0;
		cvdata->pwno[k] = 0.0;
	}

	return NR_SUCCESS;
}

// buffer call
int callnr(SEOPT opt, WAVBUFF *spdata, CAOV *cvdata)
{
	int i;	
	//static int disp;
	int enhstatus = NR_SUCCESS;

	// 잡음처리 없이 32ms delay
	if ((opt.se == 0) && (opt.it == 0))
	{
		for (i=0; i<2*BUFFSIZE; i++)
		{
			spdata->enhdata[i] = spdata->prev_bufdata[i];
			spdata->prev_bufdata[i] = spdata->bufdata[i];
		}
	}
	else
	{
		if (opt.mchginit < NUM_INITBUF)
		{
			// 잡음 추정 초기화
			//printf("Noise Estimation Initialization: %d\n", opt.mchginit);
			enhstatus = doenhinit(opt, spdata, cvdata);
			//disp = 0;
		}
		else
		{
			// 잡음 저감 함수			
			/*
			if (disp == 0)
			{
				printf("-Speech Enhancement Mode Setting \n");
				printf("SE OPT: %d\n", opt.se);
				printf("NE OPT: %d\n", opt.ne);
				printf("WT OPT: %d\n", opt.wt);
				printf("IT OPT: %d\n", opt.it);
				printf("GE OPT: %f\n", opt.ge);
				disp = 1;				
			}
			*/

			enhstatus = doenhbuf(opt, spdata, cvdata);
		}		
	}	

	return enhstatus;
}
