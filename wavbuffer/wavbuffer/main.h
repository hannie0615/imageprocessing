// Data types
struct WAVEFORM
{
	WORD wFormatTag;
	WORD nChannels;
	DWORD nSamplesPerSec;
	DWORD nAvgBytesPerSec;
	WORD nBlockAlign;
	WORD wBitsPerSample;
};

struct FMT
{
	char fmtID[4];
	DWORD fmtSIZE;
	WAVEFORM fmtFORMAT;
};

// functions


unsigned short bit2us(int *b);
int wave_read(char* filename,FMT& fmt,byte*& data,int& data_size);
int chunk_code(FILE* fp_wav,char a,char b,char c,char d);