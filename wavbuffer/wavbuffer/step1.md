## sample별로 음성 데이터(byte) 수치 출력이 가능한 경우
### C언어:
-----------------


  // ConsoleApplication2.cpp : DLL 응용 프로그램을 위해 내보낸 함수를 정의합니다.
  
  #include "stdafx.h"
  #include "main.h"
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
  #include "spenhbuff.h"
  #include "buffercnt.h"
  #include "enh_const.h"
  #include "bessel.h"
  #include "wavread.h"
  // Main routine
  // int main(int argc, char **argv)
  int WAVBUFFER_API _stdcall main(int num){
	char *infilename="sp0764_computer_sn0.wav";
	//char *outfilename="sp0764_computer_sn0_enh.wav"; 
	
	FILE *fpinfile  = NULL;
	FILE *fpoutfile = NULL;
	/*variables*/
	int i=0,j=0;
	//clock_t start, end;
	FMT fmt; /*wave file format*/
	byte* data = NULL; /*wave file data*/
	int data_size; /*data size*/
	int wave_read_status = wave_read(infilename, fmt, data, data_size);
	int num_channels = 1;
	int num_sampling_rate = 16000;
	int num_bits_per_sample = 16;
	int num_bytes_per_sample = num_bits_per_sample/8;
	int numerator = (int)pow((double)2.0,(double)num_bits_per_sample-1.0);
	
	int num_sample = data_size/(num_bytes_per_sample*num_channels);
	// Sampling Rate Check (only 16000 Hz)
	int sample_freq = num_sampling_rate;
	if ((sample_freq != DEFAULT_SAMPLING_RATE) || (num_channels != 1) || (num_bytes_per_sample !=2))
	{
		cout << "Check sampling rate of the input audio. Only Work for 16000Hz Mono 16bit Audio." << endl;
		return NR_FAILURE;
	}
	// 잡음 처리 모드
	SEOPT option = {0, 1, 0, 0, 0, 1.0}; // option for speech enhancement	
	// 버퍼 처리 변수 초기화	
	WAVBUFF spdata;
	CAOV cvdata;
	initbuff(&spdata, &cvdata);
	// for output
	byte* enhau;
	enhau = new byte[data_size];
	// initialized to zero
	for(j=0; j<data_size; j++) //data_size -> 44
	{			
		enhau[j] = 0;
	}
	 
	// data에서 32ms씩 버퍼형태로 입력된다고 가정하고, 32ms 처리하도록 코딩
	int numbuffcall = (int)(floor(double(num_sample-BUFFSIZE)/BUFFSIZE))+1; // num_sample: length of input wav file
	//int numbuffcall = 1; // 임의의 숫자 가정
	int inx;
	int initcount = 0;
	int curtime = 0;
	
	for(i=0; i < numbuffcall; i++)
	{	
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
		// spdata, cvdata 
		callnr(option, &spdata, &cvdata);
		// 끊김 없는 처리 위해서 이전 결과 출력
		// 이전 32ms 처리 결과 enhdata 출력
		for(j=0; j<2*BUFFSIZE; j++)
		{			
			enhau[inx+j] = spdata.enhdata[j];
		}
	}
	int result;
	result = enhau[num];
	return result;
}


### Python: 
-----------------

"""
call dll file
"""
from ctypes import *
import os
import numpy as np
import time
# 로컬 파일
path2 = 'C:/Users/NO. 1/Documents/WorkHome/wavbuffer/Debug/wavbuffer.dll'
input = 'sp0764_computer_sn0.wav'
dll2 = windll.LoadLibrary(path2)
v = list(range(2200,2300)) # 임의로 정한 구간 (2000,2100)
for i in v:
    print(dll2.main(i))

