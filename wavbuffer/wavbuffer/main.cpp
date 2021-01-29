#include <limits>
#include <list>
#include <fstream>
#include <iostream>
#include <sstream>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>
#include <string>

using namespace std;

#include "main.h"
#include "spenhbuff.h"
#include "buffercnt.h"
#include "enh_const.h"

// Main routine
//int main(int argc, char **argv)
int main(void)
{

	char *infilename="sp0764_computer_sn0.wav";
	char *outfilename="sp0764_computer_sn0_enh.wav";
	
	FILE *fpinfile  = NULL;
	FILE *fpoutfile = NULL;
	
	/*variables*/
	int i=0,j=0;
	//clock_t start, end;

	FMT fmt; /*wave file format*/
	byte* data = NULL; /*wave file data*/
	int data_size; /*data size*/
	int wave_read_status = wave_read(infilename, fmt, data, data_size);

	/*Parse .wav file data into mono type*/
	int num_channels = fmt.fmtFORMAT.nChannels;
	int num_sampling_rate=fmt.fmtFORMAT.nSamplesPerSec;
	int num_bits_per_sample = fmt.fmtFORMAT.wBitsPerSample;
	int num_bytes_per_sample = num_bits_per_sample/8;
	int numerator = (int)pow((double)2.0,(double)fmt.fmtFORMAT.wBitsPerSample-1.0);
	
	int num_sample = data_size/(num_bytes_per_sample*num_channels);

	// Sampling Rate Check (only 16000 Hz)
	int sample_freq = fmt.fmtFORMAT.nSamplesPerSec;
	if ((sample_freq != DEFAULT_SAMPLING_RATE) || (num_channels != 1) || (num_bytes_per_sample !=2))
	{
		cout << "Check sampling rate of the input audio. Only Work for 16000Hz Mono 16bit Audio." << endl;
		return NR_FAILURE;
	}

	printf("Input  speech: %s\n", infilename);
	printf("Output speech: %s\n", outfilename);
	printf("data_size: %d\n", data_size);
	printf("num_sample: %d\n", num_sample);
	printf("num_sampling_rate: %d\n", fmt.fmtFORMAT.nSamplesPerSec);
	printf("num_channel: %d\n", fmt.fmtFORMAT.nChannels);
	printf("num_bits_per_sample: %d\n", fmt.fmtFORMAT.wBitsPerSample);

	// 잡음 처리 모드
	SEOPT option = {0, 1, 0, 0, 0, 1.0}; // option for speech enhancement	

	printf("-Speech Enhancement Mode Setting (Initial Setting) \n");
	printf("SE OPT: %d\n", option.se);
	printf("NE OPT: %d\n", option.ne);
	printf("WT OPT: %d\n", option.wt);
	printf("IT OPT: %d\n", option.it);
	printf("GE OPT: %f\n", option.ge);

	// 버퍼 처리 변수 초기화	
	WAVBUFF spdata;
	CAOV cvdata;

	initbuff(&spdata, &cvdata);

	// for output
	byte* enhau;
	enhau = new byte[data_size];

	// initialized to zero
	for(j=0; j<data_size; j++)
	{			
		enhau[j] = 0;
	}
			
	// data에서 32ms씩 버퍼형태로 입력된다고 가정하고, 32ms 처리하도록 코딩
	int numbuffcall = (int)(floor(double(num_sample-BUFFSIZE)/BUFFSIZE))+1; // num_sample: length of input wav file
	int inx;
	int initcount = 0;
	int curtime = 0;
	for(i=0; i < numbuffcall; i++)
	{
		// 모드 변경 예시 1 (4초 시점에서 모드 변경)			
		/*
		if (i==125) 
		{
		    curtime = (int)(float(i*BUFFSIZE)/DEFAULT_SAMPLING_RATE);
			printf("\nSpeech Enhancement option changed at %d sec.\n", curtime);
			option.se = 2; 
			option.wt = 1; 
			//initcount = 0;

			printf("-Speech Enhancement Mode Setting \n");
			printf("SE OPT: %d\n", option.se);
			printf("NE OPT: %d\n", option.ne);
			printf("WT OPT: %d\n", option.wt);
			printf("IT OPT: %d\n", option.it);
			printf("GE OPT: %f\n", option.ge);
		}
		*/

		// 모드 변경 예시 2 (4초 시점에서 gain adjustment)		
		/*
		if (i==125) 
		{
			option.ge = 1.7;
			curtime = (int)(float(i*BUFFSIZE)/DEFAULT_SAMPLING_RATE);
			printf("\nSpectral gain (GE OPT) adjusted to %3.2f at %d sec.\n", option.ge, curtime);
			printf("-Speech Enhancement Mode Setting \n");
			printf("SE OPT: %d\n", option.se);
			printf("NE OPT: %d\n", option.ne);
			printf("WT OPT: %d\n", option.wt);
			printf("IT OPT: %d\n", option.it);
			printf("GE OPT: %f\n", option.ge);
		}
		*/
		
		
		// 모드 변경 후 첫 3개 버퍼에서 잡음추정
		if (initcount < NUM_INITBUF+1)
		{
			option.mchginit = initcount;
			initcount = initcount+1;
		}		

		inx = i*2*BUFFSIZE; //data는 byte형이고, 음성은 16비트 인코딩이므로 음성 1샘플은 2byte임. 따라서 (2*BUFFSIZE)만큼 읽어야 32ms임

		// 현재 32ms 버퍼 입력
		for(j=0; j<2*BUFFSIZE; j++)
		{			
			spdata.bufdata[j] = data[inx+j];
		}

		callnr(option, &spdata, &cvdata);

		// 끊김 없는 처리 위해서 이전 결과 출력
		// 이전 32ms 처리 결과 enhdata 출력
		for(j=0; j<2*BUFFSIZE; j++)
		{			
			enhau[inx+j] = spdata.enhdata[j];
		}
	}
					
	// Wav write
	// infile & outfile share the same header
	fpinfile = fopen(infilename,"rb");
	byte header[44];
	fread(header,sizeof(byte),44,fpinfile);	
	fclose(fpinfile);

	fpoutfile = fopen(outfilename,"wb");
	fwrite(header,sizeof(byte),44,fpoutfile);	
	fwrite(enhau,sizeof(byte),data_size,fpoutfile);
	//fwrite(data,2*sizeof(byte),num_sample,fpoutfile);
	fclose(fpoutfile);

	return NR_SUCCESS;
}

int wave_read(char* filename,FMT& fmt,byte*& data,int& data_size){
	/*read wav*/
	FILE* fp_wav = NULL;
	fp_wav = fopen(filename,"rb");
	if ( fp_wav == NULL ){
	//	error_mmp_ar(MMP_AR_FINGERPRINT_EXT_WAVE_READ_ERROR_FILE_OPEN);
		//return MMP_AR_FINGERPRINT_EXT_WAVE_READ_ERROR_FILE_OPEN;
	}

	/*read file and check if it is 'RIFF'*/
	if ( chunk_code(fp_wav,'R','I','F','F') != NR_SUCCESS ){
	//	error_mmp_ar(MMP_AR_FINGERPRINT_EXT_WAVE_READ_ERROR_WAVE_READ_ERROR);
		//return MMP_AR_FINGERPRINT_EXT_WAVE_READ_ERROR_WAVE_READ_ERROR;
	}

	DWORD remaining;
	/*Read remaining*/
	if ( !fread((void*)&remaining,sizeof(DWORD),1,fp_wav) ){
		fclose(fp_wav);
	//	error_mmp_ar(MMP_AR_FINGERPRINT_EXT_WAVE_READ_ERROR_WAVE_READ_ERROR);
		//return MMP_AR_FINGERPRINT_EXT_WAVE_READ_ERROR_WAVE_READ_ERROR;
	}

	/*read file and check if it is 'WAVE'*/
	if ( chunk_code(fp_wav,'W','A','V','E') != NR_SUCCESS ){
	//	error_mmp_ar(MMP_AR_FINGERPRINT_EXT_WAVE_READ_ERROR_WAVE_READ_ERROR);
		//return MMP_AR_FINGERPRINT_EXT_WAVE_READ_ERROR_WAVE_READ_ERROR;
	}
	remaining -= 4;

	/*read fmt chunk*/
	if ( !fread((void*)&fmt,sizeof(FMT),1,fp_wav) ){
		fclose(fp_wav);
	//	error_mmp_ar(MMP_AR_FINGERPRINT_EXT_WAVE_READ_ERROR_WAVE_READ_ERROR);
		//return MMP_AR_FINGERPRINT_EXT_WAVE_READ_ERROR_WAVE_READ_ERROR;
	}
	remaining -= 4;/*for 'fmt '*/
	remaining -= 4;/*for fmt.fmtSIZE*/
	remaining -= fmt.fmtSIZE;/*for WAVEFORM*/

	/*check fmt*/
	if ( fmt.fmtID[0] != 'f' || fmt.fmtID[1] != 'm' || fmt.fmtID[2] != 't' ){
		fclose(fp_wav);
	//	error_mmp_ar(MMP_AR_FINGERPRINT_EXT_WAVE_READ_ERROR_WAVE_READ_ERROR);
		//return MMP_AR_FINGERPRINT_EXT_WAVE_READ_ERROR_WAVE_READ_ERROR;
	}
	/*when there is an additional bytes in fmt chunk*/
	if (fmt.fmtSIZE > sizeof(FMT)){
		byte* tmp = new byte[fmt.fmtSIZE-sizeof(WAVEFORM)];
		if ( !fread((void*)tmp,fmt.fmtSIZE-sizeof(WAVEFORM),1,fp_wav) || tmp==NULL ){
			fclose(fp_wav);
			return NR_IO_WAVE_READ_WAV_READ_ERROR;
		}
		if(tmp) delete[] tmp;
	}

	/*check supported fmt*/
	if ( fmt.fmtFORMAT.wFormatTag != 1){/*1 is for PCM*/
		fclose(fp_wav);
		return NR_IO_WAVE_READ_WAV_READ_ERROR;
	}

	/*read file and check if it is 'DATA'*/
	if ( chunk_code(fp_wav,'d','a','t','a') != NR_SUCCESS ){
		return NR_IO_WAVE_READ_WAV_READ_ERROR;
	}
	remaining -= 4;

	/*now read data of data chunk*/
	/*check remaining*/
	if ( remaining < 0 ){
		fclose(fp_wav);
		return NR_IO_WAVE_READ_WAV_READ_ERROR;
	}
	
	/*Get lenght of sound data*/
	DWORD length;
	/*Read length*/
	if ( !fread((void*)&length,sizeof(DWORD),1,fp_wav) ){
		fclose(fp_wav);
		return NR_IO_WAVE_READ_WAV_READ_ERROR;
	}

	data = new byte[length];
	if (data == NULL){/*Insufficient memory*/	
		return NR_IO_WAVE_READ_WAV_MEMORY_ASSIGN;
	}

	data_size = fread((void*)data,sizeof(byte),length,fp_wav);

	if ( data_size < 0 ){
		data_size = 0;
		delete data;
		fclose(fp_wav);	
		return NR_IO_WAVE_READ_WAV_MEMORY_ASSIGN;
	}
	
	fclose(fp_wav);
	return NR_SUCCESS;
}

int chunk_code(FILE* fp_wav,char a,char b,char c,char d){
	char chunk_type[4];
	/*Read 'abcd'*/
	while(1){
		if ( fread((void*)chunk_type,sizeof(byte),1,fp_wav) ){
			if ( a == chunk_type[0] ){
				if ( fread((void*)chunk_type,sizeof(byte),3,fp_wav) ){
					chunk_type[3] = chunk_type[2];
					chunk_type[2] = chunk_type[1];
					chunk_type[1] = chunk_type[0];
					chunk_type[0] = a;
					break;
				}
				else{
					return NR_IO_CHUNK_CODE_FILE_READ;
				}

			}
			else{
				continue;
			}
		}
		else{
			return NR_IO_CHUNK_CODE_FILE_READ;
		}
	}

	/*Check match between read code and 'abcd'*/
	if (chunk_type[0] != a ||	chunk_type[1] != b || chunk_type[2] != c ||	chunk_type[3] != d ){
		return NR_IO_CHUNK_CODE_TYPE_ERR;
	}

	return NR_SUCCESS;
}



// Bit to unsigned short
unsigned short bit2us(int *b)
{
	int i;
	unsigned short us;

	us = 0;
	for (i = 0; i < 16; i++) us += b[i] * 0x8000 >> i; 

	return us;
}



