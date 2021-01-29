#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "enh_const.h"
#include "spenhbuff.h"
#include "spfilt.h"
#include "bessel.h"

// noise reduction using buffer
int doenhinit(SEOPT opt, WAVBUFF *spdata, CAOV *cvdata)
{
	int i, k, q;
	int lenf = FFTSIZE/2+1;
	
	double ar[FFTSIZE];	
	double ai[FFTSIZE];
	double fr[FFTSIZE];	
	double fi[FFTSIZE];

	double spec[FFTSIZE];		

	int num_sample = BUFFSIZE+BUFFSIZE/2; // 현재 32ms(512 samples)+ 이전 16ms(256 samples) for 16kHz sampling rate
	int num_frame = (int)(floor(double(num_sample-FRSIZE)/FRSTEP))+1;

	// 메모리를 잡고 샘플을 읽어 들인다.
	float* aublk = NULL;
	aublk = new float[num_sample];    
	if (aublk == NULL){/*Insufficient memory*/	
		return NR_NOISE_DOENHINIT_MEMORY_ASSIGN;
	}

	int num_bytes_per_sample = spdata->num_bits_per_sample/8;
	int numerator = (int)pow((double)2.0,(double)spdata->num_bits_per_sample-1.0);
	// 1ch로 가정
	
	// from the previous buffer
	short* short_data1 = reinterpret_cast<short*>(spdata->prev_bufdata);
	for(i=0; i<BUFFSIZE/2; i++)
	{
		aublk[i] = (float)(short_data1[i+BUFFSIZE/2]);
		aublk[i] /= (float)(numerator);
	}

	// from the current buffer
	short* short_data2 = reinterpret_cast<short*>(spdata->bufdata);	
	//for(i=0; i<BUFFSIZE; i++)
	for(i=BUFFSIZE/2; i<num_sample; i++)
	{
		aublk[i] = (float)(short_data2[i-BUFFSIZE/2]);
		//aublk[i+BUFFSIZE/2] = (float)(short_data2[i]);
		aublk[i] /= (float)(numerator);
	}	
	
	// 새로 시작하거나 option 변화가 있을 경우
	if (opt.mchginit == 0)
	{
		// 데이터 ar, ai로 복사 및 Windowing
		for(k=0; k<FRSIZE; k++)
		{
			ar[k] = (double)aublk[BUFFSIZE/2+k]*cwin[k];
			ai[k] = 0.0;

			fr[k] = 0.0;
			fi[k] = 0.0;
			spec[k] = 0.0;
		}

		// zero padding
		for(k=FRSIZE; k<FFTSIZE; k++)
		{
			ar[k] = 0.0;
			ai[k] = 0.0;

			fr[k] = 0.0;
			fi[k] = 0.0;
			spec[k] = 0.0;
		}
					
		// Forward FFT & Power Spectrum
		fft1(ar, ai, FFTSIZE, 0);
		for(q=0; q<lenf; q++)
		{
			//spec[q] = fr[q]*fr[q]+fi[q]*fi[q];
			spec[q] = ar[q]*ar[q]+ai[q]*ai[q];
			cvdata->noise_ps[q] = spec[q];
		}		
	}
	else	
	{
		for(i=0; i<num_frame; i++)
		{
			// 데이터 ar, ai로 복사 및 Windowing
			for(k=0; k<FRSIZE; k++)
			{
				ar[k] = (double)aublk[i*FRSTEP+k]*cwin[k];
				ai[k] = 0.0;

				fr[k] = 0.0;
				fi[k] = 0.0;
				spec[k] = 0.0;
			}

			// zero padding
			for(k=FRSIZE; k<FFTSIZE; k++)
			{
				ar[k] = 0.0;
				ai[k] = 0.0;

				fr[k] = 0.0;
				fi[k] = 0.0;
				spec[k] = 0.0;
			}

			// Forward FFT & Power Spectrum			
			fft1(ar, ai, FFTSIZE, 0);
			for(q=0; q<lenf; q++)
			{
				//spec[q] = fr[q]*fr[q]+fi[q]*fi[q];
				spec[q] = ar[q]*ar[q]+ai[q]*ai[q];
				cvdata->noise_ps[q] = cvdata->noise_ps[q]+spec[q];
			}
		}
	}

	if (opt.mchginit == NUM_INITBUF-1)
	{
		int numf = NUM_INITBUF*2-1; // 첫번째 버퍼는 1 프레임만 처리
		for(q=0; q<lenf; q++)
		{
			cvdata->noise_ps[q] = cvdata->noise_ps[q]/numf; // initialization for noise estimation
			cvdata->prev_specenh[q] = cvdata->noise_ps[q]; // initialization for ksi estimation
			cvdata->pwsp[q] = cvdata->noise_ps[q]; // initialization for intelligibility
			cvdata->pwno[q] = cvdata->noise_ps[q]; // initialization for intelligibility

			//printf("noise_ps %d, %f\n", q, cvdata->noise_ps[q]);
		}

		// martin_nspow initialization
		if (opt.ne == 0)
		{
			initialize_msparm(cvdata->noise_ps, &(cvdata->msnoise));
		}

		for (i=0; i<BUFFSIZE; i++)
		{
			spdata->nextenh[i] = spdata->bufdata[i];
		}
	}

	// update enhdata & prev_bufdata
	for (i=0; i<2*BUFFSIZE; i++)
	{
		spdata->enhdata[i] = spdata->prev_bufdata[i];
		spdata->prev_bufdata[i] = spdata->bufdata[i];
	}
	
	if(aublk) delete[] aublk;

	return NR_SUCCESS;
}

// noise reduction using buffer
int doenhbuf(SEOPT opt, WAVBUFF *spdata, CAOV *cvdata)
{
	int i, k, q;
	int lenf = FFTSIZE/2+1;
	
	double ar[FFTSIZE];	
	double ai[FFTSIZE];
	double fr[FFTSIZE];	
	double fi[FFTSIZE];

	double spec[FFTSIZE];	
	double angle[FFTSIZE];	
	double xold[FRSTEP];
	//double xiw[FRSTEP];

	for(k=0; k < FRSTEP; k++)
	{
		xold[k] = cvdata->xold[k];
	}

	int num_sample = BUFFSIZE+BUFFSIZE/2; // 현재 32ms(512 samples)+ 이전 16ms(256 samples) for 16kHz sampling rate
	int num_frame = (int)(floor(double(num_sample-FRSIZE)/FRSTEP))+1;

	double gammak[FFTSIZE];
	double ksi[FFTSIZE];
	double ksi_min = pow(10.0, -2.5);	
	double spegain[FFTSIZE];	
	//double prev_specenh[FFTSIZE];
	double mag[FFTSIZE];
	double ddalp = 0.98;
	double gse = opt.ge;
	double tmp, wtz, wtadd, zeta;

	// variables for intelligibility
	double itgain[FFTSIZE];
	double sii;
	int Bfc[NUM_ITBAND+1] = {9, 11, 14, 17, 22, 28, 35, 45, 57, 71, 91, 114, 142, 179, 228, 285, 359, 456, 512};

	// 메모리를 잡고 샘플을 읽어 들인다.
	float* aublk = NULL;
	aublk = new float[num_sample];    
	if (aublk == NULL){/*Insufficient memory*/	
		return NR_NOISE_DOENHBUF_MEMORY_ASSIGN;
	}
	
	float enhaublk[BUFFSIZE];
	short enhshort[BUFFSIZE];

	int num_bytes_per_sample = spdata->num_bits_per_sample/8;
	int numerator = (int)pow((double)2.0,(double)spdata->num_bits_per_sample-1.0);
	// 1ch로 가정
	
	// from the previous buffer
	short* short_data1 = reinterpret_cast<short*>(spdata->prev_bufdata);
	for(i=0; i<BUFFSIZE/2; i++)
	{
		aublk[i] = (float)(short_data1[i+BUFFSIZE/2]);
		aublk[i] /= (float)(numerator);
	}

	// from the current buffer
	short* short_data2 = reinterpret_cast<short*>(spdata->bufdata);	
	//for(i=0; i<BUFFSIZE; i++)
	for(i=BUFFSIZE/2; i<num_sample; i++)
	{
		aublk[i] = (float)(short_data2[i-BUFFSIZE/2]);
		//aublk[i+BUFFSIZE/2] = (float)(short_data2[i]);
		aublk[i] /= (float)(numerator);
	}	

	for(i=0; i < num_frame; i++)
	{
		// 데이터 ar, ai로 복사 및 Windowing
		for(k=0; k<FRSIZE; k++)
		{
			ar[k] = (double)aublk[i*FRSTEP+k]*cwin[k];
			ai[k] = 0.0;

			fr[k] = 0.0;
			fi[k] = 0.0;
			spec[k] = 0.0;
		}

		// zero padding
		for(k=FRSIZE; k<FFTSIZE; k++)
		{
			ar[k] = 0.0;
			ai[k] = 0.0;

			fr[k] = 0.0;
			fi[k] = 0.0;
			spec[k] = 0.0;
		}
					
		// Forward FFT & Power Spectrum
		fft1(ar, ai, FFTSIZE, 0);
		for(q=0; q<lenf; q++)
		{
			spec[q] = ar[q]*ar[q]+ai[q]*ai[q];
		}			

		for(q=0; q<FFTSIZE; q++)
		{
			angle[q] = atan2(ai[q],ar[q]);
			//angle[q] = atan2(fi[q],fr[q]);
		}
		
		int nest_status;
		if (opt.ne == 0)
		{
			nest_status = martin_nspow(spec, lenf, cvdata->noise_ps, &(cvdata->msnoise));
		}
		else if (opt.ne == 1)
		{
			nest_status = ummse_nspow(spec, lenf, cvdata->noise_ps, cvdata->prev_ph);
		}
		else
		{
			printf("unsupported noise estimation option\n");
			return NR_NOISE_DOENHBUF_NOISE_OPTION;
		}

		if (nest_status != NR_SUCCESS)
		{
			return nest_status;
		}
		
		// posteriori & priori SNR estimation
		for(q=0; q<lenf; q++)
		{
			tmp = spec[q]/(cvdata->noise_ps[q]+EPS);
			gammak[q] = getmin(tmp, 40); // posteriori SNR
		}
					
		if (opt.wt == 0)
		{
			for(q=0; q<lenf; q++)
			{
				tmp = (cvdata->prev_specenh[q])/(cvdata->noise_ps[q]+EPS);
				ksi[q] = ddalp*tmp+(1-ddalp)*getmax(gammak[q]-1,0);				
				ksi[q] = getmax(ksi[q], ksi_min);
			}
		}
		else if (opt.wt == 1)
		{
			for(q=0; q<lenf; q++)
			{
				zeta = (cvdata->prev_specenh[q])/(cvdata->noise_ps[q]+EPS);
				tmp = getmax(gammak[q]-1,0);
				wtz = tmp/(tmp+zeta+1);
				wtadd = getmin(1-wtz, wtz);

				if (wtadd < 0.1)
				{
					wtadd = ddalp*(wtadd/0.1);
				}
				else
				{
					wtadd = ddalp;
				}

				ksi[q] = wtadd*zeta+(1-wtadd)*tmp;
				ksi[q] = getmax(ksi[q], ksi_min);
			}				
		}
		else
		{
			printf("unsupported a priori snr estimation option\n");
			return NR_NOISE_DOENHBUF_SNR_OPTION;
		}
	
		// Spectral gain function for noise reduction
		if (opt.se == 1)
		{
			mmse_seg(ksi, gammak, gse, lenf, spegain);
		}
		else if (opt.se == 2)
		{
			wiener_seg(ksi, gse, lenf, spegain);
		}
		else if (opt.se == 0)
		{
			for(q=0; q<lenf; q++)
			{
				spegain[q] = 1.0;
			}
		}
		else
		{
			printf("unsupported Spectral gain function option\n");
			return NR_NOISE_DOENHBUF_SPGAIN_OPTION;
		}	

		// Spectral gain function for intelligibility enhancement
		if (opt.it == 0)
		{
			for(q=0; q<lenf; q++)
			{
				itgain[q] = 1.0;
			}
		}
		else if (opt.it == 1)
		{			
			compsii(spec, cvdata->noise_ps, cvdata->pwsp, cvdata->pwno, lenf, Bfc, &sii);			
			compitgain(sii, cvdata->pwsp, lenf, Bfc, ksi, gammak, itgain);
		}
		else
		{
			printf("unsupported intelligibility enhancement option\n");
			return NR_NOISE_DOENHBUF_INTG_OPTION;
		}			

		// Speech Enhancement by multiplying gain
		for(q=0; q<lenf; q++)
		{		
			mag[q] = spegain[q]*sqrt(itgain[q])*sqrt(spec[q]);

			cvdata->prev_specenh[q] = mag[q]*mag[q];
		}

		// fold magitude
		for(q=FFTSIZE/2+1; q<FFTSIZE; q++)
		{
			mag[q] = mag[FFTSIZE-q];
			spec[q] = spec[FFTSIZE-q];
		}				

		// reconstruct fr & fi from spec and angle
		for(q=0; q<FFTSIZE; q++)
		{
			//fr[q] = mag[q]*ar[q]/(sqrt(spec[q])+EPS);
			//fi[q] = mag[q]*ai[q]/(sqrt(spec[q])+EPS);
			fr[q] = mag[q]*cos(angle[q]);
			fi[q] = mag[q]*sin(angle[q]);
		}
		
		// Inverse FFT				
		fft1(fr, fi, FFTSIZE, 1);		

		// overlap ADD (Note: 50% overlap 가정함)
		for(k=0; k<FRSTEP; k++)
		{
			enhaublk[i*FRSTEP+k] = (float)(xold[k]+fr[k]);
			xold[k] = fr[FRSTEP+k];
			//enhaublk[i*FRSTEP+k] = xold[k]+ar[k];			
			//xold[k] = ar[FRSTEP+k];

			if (enhaublk[i*FRSTEP+k] > 1.0)
			{
				enhaublk[i*FRSTEP+k] = 1.0;
			}

			if (enhaublk[i*FRSTEP+k] < -1.0)
			{
				enhaublk[i*FRSTEP+k] = -1.0;
			}

			enhshort[i*FRSTEP+k] = (short)(numerator*enhaublk[i*FRSTEP+k]);
		}
	}

	// update xold
	for(k=0; k < FRSTEP; k++)
	{
		cvdata->xold[k] = xold[k];
	}
	
	// update prev_bufdata
	for (i=0; i<2*BUFFSIZE; i++)
	{		
		spdata->prev_bufdata[i] = spdata->bufdata[i];
	}

	byte* enhbyte = reinterpret_cast<byte*>(enhshort);

	// update enhdata	
	for (i=0; i<BUFFSIZE; i++)
	{
		spdata->enhdata[i] = spdata->nextenh[i];
		spdata->nextenh[i] = enhbyte[i+BUFFSIZE];
	}

	for (i=BUFFSIZE; i<2*BUFFSIZE; i++)
	{
		spdata->enhdata[i] = enhbyte[i-BUFFSIZE];
	}
		
	if(aublk) delete[] aublk;

	return NR_SUCCESS;
}


// MMSE noise reduction
int mmse_seg(double *ksi, double *gammak, double gse, int lenf, double *spegain)
{
	int q;
	double vk, bess0, bess1;
	double tmpa, tmpb, tmpc;

	// MMSE gain function
	for(q=0; q<lenf; q++)
	{
		vk = ksi[q]*gammak[q]/(gse+ksi[q]+EPS);
		//vk = ksi[q]*gammak[q]/(1+ksi[q]);

		// to prevent the overflow of Bessel functions
		if (vk > 1000) 
		{
			vk = 1000;			
			//printf("vk warning (just for checking)\n");
		}

		bess0 = bessi0(vk/2.0);
		bess1 = bessi1(vk/2.0);
		tmpc = exp(-0.5*vk);
		tmpa = (0.8862*pow(vk,0.5)*tmpc)/gammak[q];
		tmpb = (1+vk)*bess0+vk*bess1;
		spegain[q] = tmpa*tmpb;		
	}

	return NR_SUCCESS;
}


// Wiener noise reduction
int wiener_seg(double *ksi, double gse, int lenf, double *spegain)
{
	int q;

	// Wiener gain function
	for(q=0; q<lenf; q++)
	{
		spegain[q] = ksi[q]/(gse+ksi[q]+EPS);
		//spegain[q] = ksi[q]/(1+ksi[q]);
		//spegain[q] = 1.0;
	}
	
	return NR_SUCCESS;
}


// FFT
int fft1(double *ar, double *ai, int n, int flag)
{
	int i,j,k,it,xp,xp2,j1,j2,iter;
	double sign,w,wr,wi,dr1,dr2,di1,di2,tr,ti,arg;

	if (n<2) return(999);
	iter = (int)(log10((double)n)/log10(2.0));
	j=1;
	for(i=0;i<iter;i++)
		j*=2;
	if(fabs((float)n-(float)j)>1.0e-6)
		return NR_FAILURE;

	sign=((flag ==1) ? 1.0 : -1.0);
	xp2=n;
	for(it=0;it<iter;it++) 
	{
		xp=xp2;
		xp2/=2;
		w=PI/xp2;
		for(k=0;k<xp2;k++)
		{
			arg=k*w;
			wr=cos(arg);
			wi=sign*sin(arg);
			i=k-xp;
			for(j=xp;j<=n;j+=xp)  
			{
				j1=j+i;
				j2=j1+xp2;
				dr1=ar[j1];
				dr2=ar[j2];
				di1=ai[j1];
				di2=ai[j2];
				tr=dr1-dr2;
				ti=di1-di2;
				ar[j1]=dr1+dr2;
				ai[j1]=di1+di2;
				ar[j2]=tr*wr-ti*wi;
				ai[j2]=ti*wr+tr*wi;
			}
		}
	}

	j1=n/2;
	j2=n-1;
	j=1;

	for(i=1;i<=j2;i++)  
	{
		if(i<j)	
		{
			tr=ar[j-1];
			ti=ai[j-1];
			ar[j-1]=ar[i-1];
			ai[j-1]=ai[i-1];
			ar[i-1]=tr;
			ai[i-1]=ti;
		}
		k=j1;
		while(k<j)  
		{
			j-=k;
			k/=2;
		}
		j+=k;
	}
	w=n;
	if(flag == 1)	
	{
		for(i=0;i<n;i++)  
		{
			ar[i]/=w;
			ai[i]/=w;
		}
	}
	return NR_SUCCESS;
}

int initialize_msparm(double *noise_ps, MSPARM *msest)
{
	int i, j;
	double tmp=0.0;

	msest->alpha_corr = 0.96;
	msest->subwc = 2;
	msest->u = 0;

	for (i=0; i<LENE; i++)
	{
		msest->alpha[i] = 0.0;
		msest->P[i] = noise_ps[i];
		msest->Pbar[i] = noise_ps[i];
		msest->Psqbar[i] = noise_ps[i];
		msest->Pmin[i] = noise_ps[i];
		msest->actmin[i] = noise_ps[i];
		msest->actmin_sub[i] = noise_ps[i];
		msest->lminflag[i] = 0;

		if (tmp<noise_ps[i])
		{
			tmp = noise_ps[i];
		}
	}

	for (i=0; i<LENE; i++)
	{
		for (j=0; j<UM; j++)
		{
			msest->minact[i][j] = tmp;
		}
	}

	return NR_SUCCESS;
}

int martin_nspow(double *spec, int len, double *prev_noise_ps, MSPARM *msest)
{
	// fixed parameters for noise estimation
	int D = 150;
	int V = 15;
	double Av = 2.12;
	double alpha_max = 0.96;
	double alpha_min = 0.3;
	double beta_max = 0.8;
	double M_D = 0.905;
	double H_D = 4.1;
	double M_V = 0.668;
	double H_V = 1.55;

	int i, j;
	double tmpa, tmpb, tmpc, tmpd;
	double alpha_corr_t;
	double alpha, bet;
	double qeqinv, qeq, qeqtild, qeqtildsub;
	double bmin[LENE], bminsub[LENE];
	double qinvbar, Bc;
	
	// calculating the optimal smoothing correction factor
	tmpa = 0.0;
	tmpb = 0.0;
	for (i=0; i<len; i++)
	{
		tmpa = tmpa+msest->P[i];
		tmpb = tmpb+spec[i];
	}
	tmpc = tmpa/(tmpb+EPS)-1;
	alpha_corr_t = 1/(1+tmpc*tmpc);
	msest->alpha_corr = 0.7*msest->alpha_corr+0.3*getmax(alpha_corr_t, 0.7);

	//printf("alpha_corr: %lf\n", msest->alpha_corr);
	//printf("len: %d\n", len);

	// calculating the optimal smoothing factor & smoothed periodogram
	tmpd = 0.0;
	for (i=0; i<len; i++)
	{
		tmpc = (msest->P[i])/(prev_noise_ps[i]+EPS) -1;
		tmpa = (alpha_max*(msest->alpha_corr))/(1+tmpc*tmpc);
		alpha = getmax(tmpa, 0.3);
		bet = getmin(alpha*alpha, beta_max);
		
		msest->P[i] = alpha*msest->P[i]+(1-alpha)*spec[i];
		msest->Pbar[i] = bet*msest->Pbar[i]+(1-bet)*msest->P[i];
		msest->Psqbar[i] = bet*msest->Psqbar[i]+(1-bet)*(msest->P[i])*(msest->P[i]);

		tmpb = fabs(msest->Psqbar[i]-(msest->Pbar[i])*(msest->Pbar[i]));
		qeqinv = tmpb/(2.0*prev_noise_ps[i]*prev_noise_ps[i]+EPS);
		qeqinv = getmin(qeqinv, 0.5);
		qeq = 1/(qeqinv+EPS);
		tmpd = tmpd+qeqinv;
		//tmpd = tmpd+1/(qeq+EPS);
		qeqtild = (qeq-2*M_D)/(1-M_D);
		qeqtildsub = (qeq-2*M_V)/(1-M_V);
		bmin[i] = 1+(D-1)*2.0/qeqtild;
		bminsub[i] = 1+(V-1)*2.0/qeqtildsub;		
	}
	qinvbar = tmpd/len;
	Bc = 1+Av*sqrt(qinvbar);

	//printf("qinvbar: %lf\n", qinvbar);
	//printf("Bc: %lf\n", Bc);

	// calculation of actmin(i,k) and actmin_sub(i,k)
	for (i=0; i<len; i++)
	{
		tmpa = (msest->P[i])*bmin[i]*Bc;
		if (tmpa < msest->actmin[i])
		{
			msest->actmin_sub[i] = (msest->P[i])*bminsub[i]*Bc;
			msest->actmin[i] = tmpa;

			if (msest->subwc==V)
			{
				msest->lminflag[i] = 0;
			}
			else
			{
				msest->lminflag[i] = 1;
			}
		}		
	}

	// noise est.	
	int uval = msest->u;
	double noise_slope_max;
	if (msest->subwc==V)
	{		
		// calculation of Pmin_u the minimum of the last U stored values of actmin
		for (i=0; i<len; i++)
		{
			// storing the value of actmin(i,k)
			msest->minact[i][uval] = msest->actmin[i];

			tmpa = msest->minact[i][0];
			for (j=1; j<UM; j++)
			{
				if (tmpa > msest->minact[i][j])
				{
					tmpa = msest->minact[i][j];
				}
			}
			msest->Pmin[i] = tmpa;
		}

		// calculation of noise slope max
		if (qinvbar < 0.03)
		{
			noise_slope_max=8;
		}
		else if (qinvbar < 0.05)
		{
			noise_slope_max=4;
		}
		else if (qinvbar < 0.06)
		{
			noise_slope_max=2;
		}
		else
		{
			noise_slope_max=1.2;		
		}

		// update Pmin_u if the minimum falls within the search range
		for (i=0; i<len; i++)
		{
			if ( (msest->lminflag[i] > 0) && (msest->actmin_sub[i]<(noise_slope_max*msest->Pmin[i])) && (msest->actmin_sub[i]>msest->Pmin[i]) )
			{
				msest->Pmin[i] = msest->actmin_sub[i];

				for (j=0; j<UM; j++)
				{
					msest->minact[i][j] = msest->actmin_sub[i];
				}					
				msest->actmin[i] = msest->actmin_sub[i];
			}

			msest->lminflag[i] = 0;
			msest->actmin[i] = msest->P[i];
			msest->actmin_sub[i] = msest->P[i];
		}

		msest->subwc = 1;
		msest->u = uval+1;
		if (msest->u == UM)
		{
			msest->u = 0;
		}
	}
	else
	{
		if (msest->subwc > 1)
		{
			for (i=0; i<len; i++)
			{
				prev_noise_ps[i] = getmin(msest->actmin_sub[i], msest->Pmin[i]);
				msest->Pmin[i] = prev_noise_ps[i];
			}
		}

		msest->subwc = msest->subwc + 1;
	}

	return NR_SUCCESS;
}


int ummse_nspow(double *spec, int len, double *prev_noise_ps, double *prev_ph)
{
	double alphaPH1mean = 0.9;
	double alphaPSD = 0.8;

	// constants for a posteriori SPP
	double q = 0.5; // a priori probability of speech presence:
	double priorFact  = q/(1-q);
	double xiOptDb = 15; // optimal fixed a priori SNR for SPP estimation
	double xiOpt = pow(10, xiOptDb/10); //10.^(xiOptDb./10)
	double logGLRFact = log(1/(1+xiOpt));
	double GLRexp = xiOpt/(1+xiOpt);

	int i;
	double snrpost, glr, ph, estimate;
	double* phmean = NULL;
	phmean = new double[len];
	if (phmean == NULL){/*Insufficient memory*/
		return NR_NOISE_NSPOW_MEMORY_ASSIGN;
	}

	double *noise_ps;
	noise_ps = new double[len];
	if (noise_ps == NULL){/*Insufficient memory*/
		if(phmean) delete[] phmean;
		return NR_NOISE_NSPOW_MEMORY_ASSIGN;
	}

	for(i=0; i<len; i++)
	{
		snrpost = spec[i]/(prev_noise_ps[i]+EPS);
		//printf("sp %d, %lf\n", i+1, spec[i]);		
		//printf("pn %d, %lf\n", i+1, prev_noise_ps[i]);
		glr = priorFact*exp(getmin(logGLRFact + GLRexp*snrpost, 200));		
		ph = glr/(1+glr);
		//printf("ph %d, %lf\n", i+1, ph);
		phmean[i] = alphaPH1mean*prev_ph[i]+(1-alphaPH1mean)*ph;
		if (phmean[i] > 0.99)
		{
			ph = getmin(ph, 0.99);
		}

		estimate = ph*prev_noise_ps[i]+(1-ph)*spec[i];
		noise_ps[i] = alphaPSD *prev_noise_ps[i]+(1-alphaPSD)*estimate;
	}
		
	// update for next frame
	for(i=0; i<len; i++)
	{
		prev_ph[i] = phmean[i];
		prev_noise_ps[i] = noise_ps[i];
	}	

	if(phmean) delete[] phmean;
	if(noise_ps) delete[] noise_ps;

	return NR_SUCCESS;
}


int compsii(double *spec, double *noise_ps, double *pwsp, double *pwno, int len, int *Bfc, double *sii)
{	
	int i, j, q;
	double wtband[NUM_ITBAND] = {0.0083, 0.0095, 0.0150, 0.0289, 0.0440, 0.0578, 0.0653, 0.0711, 0.0818, 0.0844, 0.0882, 0.0898, 0.0868, 0.0844, 0.0771, 0.0527, 0.0364, 0.0185};
	double diffband[NUM_ITBAND] = {35.6359, 44.5449, 57.9084, 75.7264, 89.0899, 115.8168, 151.4528, 178.1797, 222.7247, 311.8146, 356.3595, 445.4494, 579.0842, 757.2639, 890.8987, 1158.1683, 1514.5278, 872.8103};
	double qv;
	double tmpA, tmpB, tmpC, tmpD;
	double eqsp, eqno, dif;
	double alp = 0.98;
	double powsp, powno;

	// smoothed spectrum
	for (i=0; i<len; i++)
	{
		pwsp[i] = alp*pwsp[i]+(1-alp)*spec[i];
		pwno[i] = alp*pwno[i]+(1-alp)*noise_ps[i];		
	}
	
	tmpD = 0.0;
	for (j=0; j<NUM_ITBAND; j++)
	{
		tmpA = 0.0;
		tmpB = 0.0;
		for (q=Bfc[j]; q<Bfc[j+1]; q++)
		{
			tmpA = tmpA+pwsp[q];
			tmpB = tmpB+pwno[q];
			//tmpA = tmpA+spec[q];
			//tmpB = tmpB+noise_ps[q];
		}
		tmpC = Bfc[j+1]-Bfc[j];
		tmpA = tmpA/tmpC;
		tmpB = tmpB/tmpC;

		powsp = tmpA;
		powno = tmpB;

		eqsp = 10*log10(powsp/diffband[j]);
		eqno = 10*log10(powno/diffband[j]);
		dif = (eqsp-eqno+15)/30;
		qv = getmax(getmin(dif, 1.0), 0.0);

		tmpD = tmpD+wtband[j]*qv;
	}

	*sii = tmpD;

	return NR_SUCCESS;
}

double compitgain(double sii, double *pwsp, int len, int *Bfc, double *ksi, double *gammak, double *itgain)
{
	int i;
	double tmp, logsk, vdec, vadcp, siiva;

	// sii adjustment based on voice activity detection
	tmp = 0.0;
	for (i=0; i<len; i++)
	{
		logsk = gammak[i]*ksi[i]/(1+ksi[i])-log(1+ksi[i]);
		if ((i==0) || (i==(len-1)))
		{
			tmp = tmp+logsk/2; // DC and Phi of FFT
		}
		else
		{
			tmp = tmp+logsk;
		}
	}
	vdec = 2*tmp/FFTSIZE;

	if (vdec < 0)
	{
		vadcp = 1.0;
	}
	else if (vdec < 1)
	{
		vadcp = 1-vdec;
	}
	else
	{
		vadcp = 0;
	}

	//printf("%lf\n", vdec);
	siiva = getmin(sii+vadcp, 1.0);

	// computing intelligibility gain
	int q;
	double sumpwsp, sumpsii;
	for (q=0; q<len; q++)
	{
		itgain[q] = 1.0;
	}

	sumpwsp = 0.0;
	sumpsii = 0.0;
	for (q=Bfc[0]; q<Bfc[NUM_ITBAND]; q++)
	{
		sumpwsp = sumpwsp+pwsp[q];
		sumpsii = sumpsii+pow(pwsp[q], siiva);
	}

	for (q=Bfc[0]; q<Bfc[NUM_ITBAND]; q++)
	{		
		itgain[q] = sqrt((pow(pwsp[q], siiva)/(sumpsii+EPS))*(sumpwsp/(pwsp[q]+EPS)));
	}

	return NR_SUCCESS;
}

double getmin(double a, double b)
{
	if (a<b)
	{
		return a;
	}
	else
	{
		return b;
	}
}

double getmax(double a, double b)
{
	if (a>b)
	{
		return a;
	}
	else
	{
		return b;
	}
}
