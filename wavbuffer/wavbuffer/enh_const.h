#ifndef __gen_cons 
#define __gen_cons

typedef unsigned char byte;

// parameters
#define DEFAULT_SAMPLING_RATE 16000
#define FRSIZE 512 //512->1024
#define FFTSIZE 1024  // FFT size = 2*FRSIZE 1024->2048
#define FRSTEP 256  // 50% overlap 256->512
#define NUM_NOISEFR 5  // Initial noise estimation
#define NUM_ITBAND 18  // number of band for intelligibility
#define EPS 2.22e-16
#define PI 3.14159265358979323846
#define BUFFSIZE 512 // 32ms for 16kHz sampling rate -> 16ms로 수정
#define NUM_INITBUF 3 // number of buffers for initialization
#define LENE 513
#define UM 10

#endif


// Error Code (4 bytes)
// [1byte(system)][1byte(module)][1byte(function)][1byte(error)]
// system: 0x0A (가정)
// module: 0x01 (Data I/O), 0x02 (Noise & SNR Estimation), 0x03 (Noise Reduction), 0x04 (Intelligibility Enhancement)


// Default
#define NR_SUCCESS							0x0A000001
#define NR_FAILURE							0x0A000002

// Modules
#define NR_IO								0x0A010000
#define NR_NOISE							0x0A020000

// NR_IO functions
#define NR_IO_WAVE_READ						NR_IO | 0x0A000100
#define NR_IO_CHUNK_CODE					NR_IO | 0x0A000200

// NR_NOISE functions
#define NR_NOISE_DOENHINIT					NR_NOISE | 0x0A000100
#define NR_NOISE_DOENHBUF					NR_NOISE | 0x0A000200
#define NR_NOISE_NSPOW						NR_NOISE | 0x0A000300

// Errors in NR_IO
#define NR_IO_WAVE_READ_OPEN_FILE			NR_IO_WAVE_READ  | 0x0A000001
#define NR_IO_WAVE_READ_WAV_READ_ERROR		NR_IO_WAVE_READ  | 0x0A000002
#define NR_IO_WAVE_READ_WAV_MEMORY_ASSIGN	NR_IO_WAVE_READ  | 0x0A000003
#define NR_IO_CHUNK_CODE_FILE_READ			NR_IO_CHUNK_CODE | 0x0A000004
#define NR_IO_CHUNK_CODE_TYPE_ERR			NR_IO_CHUNK_CODE | 0x0A000005

// Errors in NR_NOISE
#define NR_NOISE_DOENHINIT_MEMORY_ASSIGN	NR_NOISE_DOENHINIT  |  0x0A000001
#define NR_NOISE_DOENHBUF_MEMORY_ASSIGN		NR_NOISE_DOENHBUF   | 0x0A000002
#define NR_NOISE_DOENHBUF_NOISE_OPTION		NR_NOISE_DOENHBUF   | 0x0A000003
#define NR_NOISE_DOENHBUF_SNR_OPTION		NR_NOISE_DOENHBUF   | 0x0A000004
#define NR_NOISE_DOENHBUF_SPGAIN_OPTION		NR_NOISE_DOENHBUF   | 0x0A000005
#define NR_NOISE_DOENHBUF_INTG_OPTION		NR_NOISE_DOENHBUF   | 0x0A000006
#define NR_NOISE_NSPOW_MEMORY_ASSIGN		NR_NOISE_NSPOW		| 0x0A000007






