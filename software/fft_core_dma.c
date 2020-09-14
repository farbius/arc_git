/*
 * fft_core_dma.c
 *
 *  Created on: 25 èþë. 2020 ã.
 *      Author: dim
 */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <strings.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <sys/mman.h>

#include "arc-dma-driver.h"


/*  Init transfer SG AXI DMA S2MM function
 * 	N_packets - amount of desired packets (length of packet by default - 1 Mb)
 * 	p_len     - length of packet (1 Mb by default)
 * 	buffer    - pointer on buffer with data (PL to PS)
 * 	return -
 * 			0 - succeed transfer
 *		   -1 - error open DMA
 *		   -2 - error malloc
 */
int fft_0_dma_init(int N_packets)
{
	printf(">>dma 1\n" );

	/* Open AXI DMA device */
	printf(">>open 0\n" );
	fd_fft_0 = open("/dev/arc-dma-0v0", O_RDWR);
	if (fd_fft_0 < 0) {
		printf("Can't open device file: %s\n", "/dev/arc-dma-0v0" );
		return -1;
	}
	printf(">>open 1\n" );

	fft_0_core_buffer = (unsigned int*) mmap ( NULL, FFT_CORE_DMA_BUFFER_SIZE , PROT_READ | PROT_WRITE , MAP_SHARED, fd_fft_0, 0);
	return 0;
}


/*  Start transfer SG AXI DMA S2MM function
 * 	N_packets - amount of desired packets (length of packet by default - 1 Mb)
 */
void fft_0_dma_start(int N_packets)
{

	ioctl ( fd_fft_0, DRIVER_IOCTL_WRITE_STATE,  1 );
	ioctl ( fd_fft_0, DRIVER_IOCTL_START_STREAM, 1 );

}


/*  Start transfer SG AXI DMA S2MM function
 * 	N_packets - amount of desired packets (length of packet by default - 1 Mb)
 */
int fft_0_dma_stop( )
{

	unsigned int actual_bytes = 0;
	unsigned int  *reg = (unsigned int *)calloc ( 1, sizeof ( unsigned int ) );
	ioctl ( fd_fft_0, DRIVER_IOCTL_READ_STATE, reg );

	int time_out = 0;
	while(*reg == 0)
	{
		ioctl ( fd_fft_0, DRIVER_IOCTL_READ_STATE, reg );
		usleep(1);
		time_out ++;
		if(time_out == 99999)
		{
			printf("TimeOut is reached \n\r");
			free(reg);
			return -1; // timeout reached
		}

	}

	actual_bytes = ioctl ( fd_fft_0, DRIVER_IOCTL_READ_STATE, reg );

	free(reg);
	return actual_bytes;
}


/* Transfer AXI DMA function
 * 	N_packets - amount of desired packets (length of packet by default - 1 Mb)
 * 	buffer    - pointer on buffer with data (PL to PS)
 * 	return -
 * 			0 - succeed transfer
 *		   -1 - error   timeout
 */

void fft_0_dma_exit()
{

	munmap(fft_0_core_buffer,FFT_CORE_DMA_BUFFER_SIZE); // free adc_raw_buffer
	close(fd_fft_0);

}

////////////////////////////////////////////////////////////////////////////////////


/*  Init transfer SG AXI DMA S2MM function
 * 	N_packets - amount of desired packets (length of packet by default - 1 Mb)
 * 	p_len     - length of packet (1 Mb by default)
 * 	buffer    - pointer on buffer with data (PL to PS)
 * 	return -
 * 			0 - succeed transfer
 *		   -1 - error open DMA
 *		   -2 - error malloc
 */
int fft_1_dma_init(int N_packets)
{
	printf(">>dma 1\n" );

	/* Open AXI DMA device */
	printf(">>open 0\n" );
	fd_fft_1 = open("/dev/arc-dma-1v0", O_RDWR);
	if (fd_fft_1 < 0) {
		printf("Can't open device file: %s\n", "/dev/arc-dma-1v0" );
		return -1;
	}
	printf(">>open 1\n" );

	fft_1_core_buffer = (unsigned int*) mmap ( NULL, FFT_CORE_DMA_BUFFER_SIZE , PROT_READ | PROT_WRITE , MAP_SHARED, fd_fft_1, 0);
	return 0;
}


/*  Start transfer SG AXI DMA S2MM function
 * 	N_packets - amount of desired packets (length of packet by default - 1 Mb)
 */
void fft_1_dma_start(int N_packets)
{
	ioctl ( fd_fft_1, DRIVER_IOCTL_WRITE_STATE,  1 );
	ioctl ( fd_fft_1, DRIVER_IOCTL_START_STREAM, 1 );

}


/*  Start transfer SG AXI DMA S2MM function
 * 	N_packets - amount of desired packets (length of packet by default - 1 Mb)
 */
int fft_1_dma_stop()
{

	unsigned int actual_bytes = 0;
	unsigned int  *reg = (unsigned int *)calloc ( 1, sizeof ( unsigned int ) );
	ioctl ( fd_fft_1, DRIVER_IOCTL_READ_STATE, reg );

	int time_out = 0;
	while(*reg == 0)
	{
		ioctl ( fd_fft_1, DRIVER_IOCTL_READ_STATE, reg );
		usleep(1);
		time_out ++;
		if(time_out == 99999)
		{
			printf("TimeOut is reached \n\r");
			free(reg);
			return -1; // timeout reached
		}

	}

	actual_bytes = ioctl ( fd_fft_1, DRIVER_IOCTL_READ_STATE, reg );
	free(reg);
	return actual_bytes;

}


/* Transfer AXI DMA function
 * 	N_packets - amount of desired packets (length of packet by default - 1 Mb)
 * 	buffer    - pointer on buffer with data (PL to PS)
 * 	return -
 * 			0 - succeed transfer
 *		   -1 - error   timeout
 */

void fft_1_dma_exit()
{

	munmap(fft_1_core_buffer,FFT_CORE_DMA_BUFFER_SIZE); // free adc_raw_buffer
	close(fd_fft_1);


}



