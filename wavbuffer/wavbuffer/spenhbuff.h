#include "buffercnt.h"

int doenhinit(SEOPT opt, WAVBUFF *spdata, CAOV *cvdata);
///////////////////////////////////////////////////////////////////////
// 잡음 저감 및 명료도 개선을 위해서 3개 버퍼를 잡음 구간으로 가정하고 잡음 추정
// 
// opt [in] : 잡음 저감 옵션
// spdata [out] : 버퍼 입력 및 결과 출력을 위한 WAVBUFF 구조체
// cvdata [out] : 다음 버퍼 처리를 위한 변수들
///////////////////////////////////////////////////////////////////////

int doenhbuf(SEOPT opt, WAVBUFF *spdata, CAOV *cvdata);
///////////////////////////////////////////////////////////////////////
// 잡음 저감 및 명료도 개선 버퍼 단위로 수행
// 
// opt [in] : 잡음 저감 옵션
// spdata [in-out] : 버퍼 입력 및 결과 출력을 위한 WAVBUFF 구조체
// cvdata [in-out] : 다음 버퍼 처리를 위한 변수들
///////////////////////////////////////////////////////////////////////

int doenh(float *aublk, int num_sample, float *enhaublk, SEOPT opt);
///////////////////////////////////////////////////////////////////////
//
// 음악신호 aublk (길이 num_sample)의 잡음을 제거하여 enhaublk 출력
// 
///////////////////////////////////////////////////////////////////////

int mmse_seg(double *ksi, double *gammak, double gse, int lenf, double *spegain);
///////////////////////////////////////////////////////////////////////
// MMSE noise reduction
// 추정된 사전 SNR인 ksi와 사후 SNR인 gammak를 이용하여
// MMSE 방법으로 스펙트럼 이득 spegain 계산함
// spegain값을 잡음 음성 스펙트럼에 곱하여 잡음을 저감함
// 추가적으로 gse를 사용하여 잡음 저감 강도 조절 가능함 
// (0.5에서 2사이, 1일 경우 원래 MMSE방법임)
// 
// ksi [in] : 추정된 사전 SNR
// gammak [in] : 계산된 사후 SNR
// gse [in] : 잡음 저감 강도 (1이 기본 모드임)
// lenf [in] : 스펙트럼의 길이
// spegain [out] : MMSE 방법으로 구한 잡음 저감을 위한 스펙트럼 이득
///////////////////////////////////////////////////////////////////////

int wiener_seg(double *ksi, double gse, int lenf, double *spegain);
///////////////////////////////////////////////////////////////////////
// Wiener noise reduction
// 추정된 사전 SNR인 ksi을 이용하여 
// Wiener필터 구현하여 스펙트럼 이득 spegain 계산함
// spegain값을 잡음 음성 스펙트럼에 곱하여 잡음을 저감함
// 추가적으로 gse를 사용하여 잡음 저감 강도 조절 가능함 
// (0.5에서 2사이, 1일 경우 원래 MMSE방법임)
// 
// ksi [in] : 추정된 사전 SNR
// gse [in] : 잡음 저감 강도 (1이 기본 모드임)
// lenf [in] : 스펙트럼의 길이
// spegain [out] : Wiener필터로 구한 잡음 저감을 위한 스펙트럼 이득
///////////////////////////////////////////////////////////////////////

int fft1(double *ar, double *ai, int n, int flag);
int martin_nspow(double *spec, int len, double *prev_noise_ps, MSPARM *msest);
///////////////////////////////////////////////////////////////////////
// 최소 통계량 기반 잡음 추정 함수
// msest에 저장된 기존 최소 통계량 정보와
// 현재 프레임의 스펙트럼 spec값을 비교하여 잡음 추정 수행
// 기존 잡음 추정값인 prev_noise_ps를 새로운 잡음 추정치로 업데이트함
// 
// spec [in] : 신호 스펙트럼
// len [in] : 스펙트럼의 길이
// prev_noise_ps [in-out] : 잡음 추정치
// msest [in-out] : 잡음 추정을 위한 변수
///////////////////////////////////////////////////////////////////////

int initialize_msparm(double *noise_ps, MSPARM *msest);
///////////////////////////////////////////////////////////////////////
// 최소 통계량 기반 잡음 추정을 위한 parameter 초기화
// 
// noise_ps [in] : 초기 잡음 추정값
// msest [out] : 잡음 추정을 위한 변수
///////////////////////////////////////////////////////////////////////

int ummse_nspow(double *spec, int len, double *prev_noise_ps, double *prev_ph);
///////////////////////////////////////////////////////////////////////
// 비편향 잡음 추정 함수
// 기존 잡음 추정값인 prev_noise_ps와 잡음 확률prev_ph를 
// 현재 프레임의 스펙트럼 spec값과 비교하여 잡음 추정 수행
// 기존 잡음 추정값인 prev_noise_ps를 새로운 잡음 추정치로 업데이트함
// 
// spec [in] : 신호 스펙트럼
// len [in] : 스펙트럼의 길이
// prev_noise_ps [in-out] : 잡음 추정치
// prev_ph [in-out] : 잡음 확률
///////////////////////////////////////////////////////////////////////

int compsii(double *spec, double *noise_ps, double *pwsp, double *pwno, int len, int *Bfc, double *sii);
///////////////////////////////////////////////////////////////////////
// Speech Intelligibility Index (SII) 계산
// 스펙트럼을 ANSI 표준에 따라 18개의 밴드로 나누고,
// 현재 프레임의 스펙트럼 spec값과 잡음 추정값 noise_ps로부터 밴드 SNR 구함
// ANSI 표준의 밴드별 가중치를 곱하여 합산하여 SII 계산함
// 
// spec [in] : 신호 스펙트럼
// noise_ps [in] : 잡음 스펙트럼
// pw_sp [in] : 시간축 방향으로 smoothed 신호 스펙트럼
// pw_no [in] : 시간축 방향으로 smoothed 잡음 스펙트럼
// len [in] : 스펙트럼의 길이
// Bfc [in] : ANSI 표준 스펙트럼 밴드 경계
// sii [out] : 프레임 명료도 척도 SII값
///////////////////////////////////////////////////////////////////////

double compitgain(double sii, double *pwsp, int len, int *Bfc, double *ksi, double *gammak, double *itgain);
///////////////////////////////////////////////////////////////////////
// SII값을 바탕으로 주파수 방향으로 신호 스펙트럼 재분배 수행
// 
// sii [in] : 프레임 명료도 척도 SII값
// pw_sp [in] : 시간축 방향으로 smoothed 신호 스펙트럼
// len [in] : 스펙트럼의 길이
// Bfc [in] : ANSI 표준 스펙트럼 밴드 경계
// ksi [in] : 추정된 사전 SNR
// gammak [in] : 계산된 사후 SNR
// itgain [out] : 명료도 개선을 위한 스펙트럼 이득 
///////////////////////////////////////////////////////////////////////


double getmin(double a, double b);
double getmax(double a, double b);
