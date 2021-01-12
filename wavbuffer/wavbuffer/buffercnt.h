#ifndef __gen_buf
#define __gen_buf

#include "enh_const.h"

// option for speech enhancement
typedef struct seoption
{
	int mchginit; // mode(option) change buffer count: 잡음재추정 (값이 2 이하일 경우) 모드 변경 후 3개 버퍼에서 잡음 추정
	int se; // speech enhancement: NONE(0), MMSE(1) or Wiener(2)
	int ne; // Noise estimation  : Martin(0), Unbiased(1)
	int wt; // Priori SNR est.   : DD (0) or ADD (1)
	int it; // Intelligibility   : IT off (0) or IT on (1)
	// gain adjustment 
    // typically from 0.5 (weak noise reduction) to 2 (strong noise reduction)
    // 1: unadjusted value
	double ge; 
} SEOPT;

// parameters for buffer control
typedef struct bufcnt
{
	// parameters of input speech wavform
	int num_sampling_rate;
	int num_bits_per_sample;

	// buffer data (최대 16 bit 양자화 가정)
	int bufsize; // LENBUF (샘플 사이즈)
	unsigned char bufdata[2*BUFFSIZE]; // 현재 buffer data (입력) <2*BUFFSIZE는 소리 1샘플당 2바이트(16비트)씩 이므로>
	unsigned char prev_bufdata[2*BUFFSIZE]; // 이전 buffer data	<2*BUFFSIZE는 소리 1샘플당 2바이트(16비트)씩 이므로>
	unsigned char enhdata[2*BUFFSIZE]; // 이전 buffer data 처리 결과 <2*BUFFSIZE는 소리 1샘플당 2바이트(16비트)씩 이므로>
	unsigned char nextenh[BUFFSIZE]; // 현재 buffer data 처리 결과의 절반(16ms)은 다음에 출력
} WAVBUFF;

// parameters for minimum-statistics noise estimation
typedef struct msnoise
{
	double alpha_corr;
	int subwc;
	int u;
	double alpha[LENE];
	double P[LENE];
	double Pbar[LENE];
	double Psqbar[LENE];
	double Pmin[LENE];
	double actmin[LENE];
	double actmin_sub[LENE];
	int lminflag[LENE];
	double minact[LENE][UM];
} MSPARM;

// variables carried over for processing next buffer
typedef struct carryover
{
	MSPARM msnoise; // for MS noise estimation(최소 통계량 기반 잡음 추정을 위한 변수)		
	double noise_ps[FFTSIZE];// 잡음 추정 결과 	
	double prev_ph[FFTSIZE]; // for unbiased noise estimation(비편향 잡음 추정을 위한 변수)
	double prev_specenh[FFTSIZE]; // for ksi estimation(사전 SNR 추정을 위한 변수)
	double pwsp[FFTSIZE]; // for intelligibility enhancement(명료도 개선을 위한 변수)
	double pwno[FFTSIZE]; // for intelligibility enhancement(명료도 개선을 위한 변수)
	double xold[FRSTEP]; // for overlap-add resynthesis from FFT(overlap-add합성 위한 변수)
} CAOV;

int initbuff(WAVBUFF *spdata, CAOV *cvdata);
///////////////////////////////////////////////////////////////////////
// 버퍼 초기화
// 실제 잡음 저감 시스템 구동 시에 1번만 호출됨
// 
// spdata [out] : WAVBUFF 구조체 초기화
// cvdata[out] : CAOV 구조체 초기화
///////////////////////////////////////////////////////////////////////

int callnr(SEOPT opt, WAVBUFF *spdata, CAOV *cvdata);
///////////////////////////////////////////////////////////////////////
// 버퍼 처리를 위한 상위 함수
// doenhinit()함수를 호출하여 잡음 추정 초기화
// doenhbuf()함수를 호출하여 잡음저감 및 명료도 개선 실행
// 
// opt [in] : 잡음 저감 옵션
// spdata [in-out] : 버퍼 입력 및 결과 출력을 위한 WAVBUFF 구조체
// cvdata [in-out] : 다음 버퍼 처리를 위한 변수들
///////////////////////////////////////////////////////////////////////


#endif
