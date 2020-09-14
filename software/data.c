#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <strings.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <sys/mman.h>

#include "stdafx.h"
#include "Soft_Part_I.h"
#include "Soft_Part_II.h"
#include "Soft_Part_III.h"
#include "recognition.h"

#include <string.h>

#include "arc-dma-driver.h"

// #define BSZ	0x01400000

#define NSP	8192
// #define WITHOUT_DMA


#define BIT_STOP_TLAST  (1<<4)
#define BIT_FFT_RST		(0x8)
#define BIT_FFT_RST		(0x8)
#define BIT_FFT_WIN 	(0x2)
#define BIT_FFT_ENAPW 	(0x1)

extern st_settings 	s;
extern st_fpgadata 	*fpgadata;
extern st_Results 		*Results;
extern st_trace 		**traces;
//extern st_TracksInfo 	*TracksInfo;
extern uint 			*TracksLengths;

static char *devName = "/dev/arc-dma-1v0";


void init_fft(void);
static int fd;
unsigned int *buf;

static void mmap_set(unsigned int* virtual_address, int offset, unsigned int value) {
	virtual_address[offset>>2] = value;
}

static unsigned int mm_get(unsigned int* virtual_address, int offset) {
    return virtual_address[offset>>2];
}

int get_Nstep (char *datafile, unsigned int *dma_buffer) // функция  для работы с драйвером ДМА, ПОТОК ЗАПУЩЕН!!!
{
	 	 printf("Start FFT_0_CORE block \n\r");
	 	 int i;
	 	 int Nstep = 0;
	 	 int error = 0, fp;


#ifdef WITHOUT_DMA


		fft_1_core_buffer = (unsigned int *) malloc ( FFT_CORE_DMA_BUFFER_SIZE );
		if(fft_1_core_buffer == NULL)
				printf("Malloc DMA buffer Error \n\r");


		int dh_0v = open("/dev/mem", O_RDWR | O_SYNC);
		unsigned int* fft_0_core = mmap(NULL, 65535, PROT_READ | PROT_WRITE, MAP_SHARED, dh_0v, 0xA0020000);

		for (i = 0 ; (i < (FFT_CORE_DMA_BUFFER_SIZE/4)) ; i++ ) //
			fft_1_core_buffer[i] = 0;

		mmap_set(fft_0_core, 0x0000, BIT_FFT_RST);// reset core 0 hw , reset by 1
		mmap_set(fft_0_core, 0x0000, 0);// reset core 0 hw , reset by 1


		mmap_set(fft_0_core, 0x0010,  0); // bin 1000  )// + 15 optional (FFFF)
		mmap_set(fft_0_core, 0x0018,  NSP+1); // NSP for fft_core
		mmap_set(fft_0_core, 0x0014, 5000);// test amplitude
		mmap_set(fft_0_core, 0x0008, 0x0008); //power trs 8= 7 bit


		mmap_set(fft_0_core, 0x0000, (BIT_FFT_WIN | BIT_STOP_TLAST) ); // 1:enable FFT window, 4  -  1: stop mode for tlast



		error = fft_1_dma_init( NSP ) ; // 34 ms
							if(error < 0)
								printf("Get fft_0_dma_init Error is %d \n\r", error);


		fft_1_dma_start(NSP); // 22 ms

		mmap_set(fft_0_core, 0x0000, (BIT_FFT_WIN | BIT_STOP_TLAST | BIT_FFT_ENAPW)); // after init only this to next data transfer

		// write(valuefd,"1", 2); // 2 s
		actual_bytes = fft_1_dma_stop();
		if(actual_bytes < 0)
			   printf("Get fft_0_dma_stop Error is %d \n\r", error);
	   //  write(valuefd,"0", 2);


		printf("actual bytes is %d \n", actual_bytes);


		mmap_set(fft_0_core, 0x0000, (BIT_FFT_WIN | BIT_STOP_TLAST) ); // 1:enable FFT window, 4  -  1: stop mode for tlast


		if(actual_bytes <= 0)
			printf("actual bytes is %d", actual_bytes);
				//	exit(0);

		if((fp = fopen(datafile,"wb")) == NULL)
			printf("Unable to open data file '%s'\n",datafile);
		else {
			fwrite(fft_1_core_buffer,actual_bytes,1,fp);
			fclose(fp);
		}

		fft_1_dma_exit();

#else

		actual_bytes = 15838672;

		fft_1_core_buffer = malloc(FFT_CORE_DMA_BUFFER_SIZE);
		if(fft_1_core_buffer <= 0)
			printf("Malloc Error \n");

		if((fp = fopen(datafile,"rb")) == NULL)
					printf("Unable to open data file '%s'\n",datafile);
				else {
					fread(fft_1_core_buffer,1 ,actual_bytes,fp);
					fclose(fp);
				}
#endif



		for (i = 0 ; (i < (actual_bytes/4)) ; i++ ) //
				if ((fft_1_core_buffer[i] & 0x3FFF) == 0x0000 && (fft_1_core_buffer[i] != 0))
				{
					Nstep++;
				}

		printf("Nstep for fft_1_dma_exit is %d \n\r", Nstep);









	return Nstep;
}

int get_data (uint Nstep, char *datafile, st_fpgadata *fpgadata, unsigned int *dma_buffer) {
	int idx;
	uint i, n;
	uint dec = 0;
	uint val;


	for(idx = 0;(fft_1_core_buffer[idx] != 0) && (idx < (actual_bytes/4)); idx++) {
		// printf("fft_1_core_buffer[%d] = %d \n",idx,fft_1_core_buffer[idx] >> 14);
		switch (fft_1_core_buffer[idx] & 0x3FFF) {

			case 0x0000: {
				if(dec == 0)
					dec = ((fft_1_core_buffer[idx] >> 14) & 0x3FFFF) - 1;
				i = ((fft_1_core_buffer[idx] >> 14) & 0x3FFFF) - dec;
				}
			break;

			default:
				n = (fft_1_core_buffer[idx] & 0x3FFF);
				// printf("fft_1_core_buffer[%d] = %d \n",n,fft_1_core_buffer[idx] >> 14);
				fpgadata[i-1].Sf[n-1] 	= (float) (fft_1_core_buffer[idx] >> 14);
				if((float) (fft_1_core_buffer[idx] >> 14) > 50)
				{
				fpgadata[i-1].Begs[n-1]	= n;
				fpgadata[i-1].Fins[n-1]	= n;
				}
			break;

			case 0x3FFF:
				fpgadata[i-1].Sf[n] 	= 0;
				fpgadata[i-1].Begs[n]	= 0;
				fpgadata[i-1].Fins[n]	= 0;
			break;

		}

	}
	// close everything here

	// the part of FPGA for a bin concatenations (in MATLAB this code start from "dI = diff(des)...") 	ToDo remove later
	int k;
	uint Begs, Fins;

	for (i = 0; i < Nstep; i++) {

		k = -1;
		Begs = fpgadata[i].Begs[0];
		for (n = 1; n < sizeof(fpgadata[i].Begs)/sizeof(uint); n++) {

			if ((fpgadata[i].Begs[n] - Begs) == 1) {
				Begs = fpgadata[i].Begs[n];
				Fins = fpgadata[i].Fins[n];
				fpgadata[i].Begs[n] = 0;
				fpgadata[i].Fins[n] = 0;
				if (k == -1) k = n-1;
			} else {
				if (k != -1) {
					fpgadata[i].Fins[k] = Fins;
					k = -1;
				}
				Begs = fpgadata[i].Begs[n];
			}

		}

	}

	// shift bins to a head of array:
	for (i = 0; i < Nstep; i++) {

		k = 0;
		for (n = 0; n < sizeof(fpgadata[i].Begs)/sizeof(uint); n++) {

			if (fpgadata[i].Begs[n] != 0) {
				if (n != k) {
					fpgadata[i].Begs[k] = fpgadata[i].Begs[n];
					fpgadata[i].Fins[k] = fpgadata[i].Fins[n];
					fpgadata[i].Begs[n] = 0;
					fpgadata[i].Fins[n] = 0;

				}
				k++;
			}

		}

	}


	return Nstep;

}
