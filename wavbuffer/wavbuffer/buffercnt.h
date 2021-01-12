#ifndef __gen_buf
#define __gen_buf

#include "enh_const.h"

// option for speech enhancement
typedef struct seoption
{
	int mchginit; // mode(option) change buffer count: ���������� (���� 2 ������ ���) ��� ���� �� 3�� ���ۿ��� ���� ����
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

	// buffer data (�ִ� 16 bit ����ȭ ����)
	int bufsize; // LENBUF (���� ������)
	unsigned char bufdata[2*BUFFSIZE]; // ���� buffer data (�Է�) <2*BUFFSIZE�� �Ҹ� 1���ô� 2����Ʈ(16��Ʈ)�� �̹Ƿ�>
	unsigned char prev_bufdata[2*BUFFSIZE]; // ���� buffer data	<2*BUFFSIZE�� �Ҹ� 1���ô� 2����Ʈ(16��Ʈ)�� �̹Ƿ�>
	unsigned char enhdata[2*BUFFSIZE]; // ���� buffer data ó�� ��� <2*BUFFSIZE�� �Ҹ� 1���ô� 2����Ʈ(16��Ʈ)�� �̹Ƿ�>
	unsigned char nextenh[BUFFSIZE]; // ���� buffer data ó�� ����� ����(16ms)�� ������ ���
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
	MSPARM msnoise; // for MS noise estimation(�ּ� ��跮 ��� ���� ������ ���� ����)		
	double noise_ps[FFTSIZE];// ���� ���� ��� 	
	double prev_ph[FFTSIZE]; // for unbiased noise estimation(������ ���� ������ ���� ����)
	double prev_specenh[FFTSIZE]; // for ksi estimation(���� SNR ������ ���� ����)
	double pwsp[FFTSIZE]; // for intelligibility enhancement(��ᵵ ������ ���� ����)
	double pwno[FFTSIZE]; // for intelligibility enhancement(��ᵵ ������ ���� ����)
	double xold[FRSTEP]; // for overlap-add resynthesis from FFT(overlap-add�ռ� ���� ����)
} CAOV;

int initbuff(WAVBUFF *spdata, CAOV *cvdata);
///////////////////////////////////////////////////////////////////////
// ���� �ʱ�ȭ
// ���� ���� ���� �ý��� ���� �ÿ� 1���� ȣ���
// 
// spdata [out] : WAVBUFF ����ü �ʱ�ȭ
// cvdata[out] : CAOV ����ü �ʱ�ȭ
///////////////////////////////////////////////////////////////////////

int callnr(SEOPT opt, WAVBUFF *spdata, CAOV *cvdata);
///////////////////////////////////////////////////////////////////////
// ���� ó���� ���� ���� �Լ�
// doenhinit()�Լ��� ȣ���Ͽ� ���� ���� �ʱ�ȭ
// doenhbuf()�Լ��� ȣ���Ͽ� �������� �� ��ᵵ ���� ����
// 
// opt [in] : ���� ���� �ɼ�
// spdata [in-out] : ���� �Է� �� ��� ����� ���� WAVBUFF ����ü
// cvdata [in-out] : ���� ���� ó���� ���� ������
///////////////////////////////////////////////////////////////////////


#endif
