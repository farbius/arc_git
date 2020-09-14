/*
 * arc-dma-driver.h
 *
 *  Created on: 25 èþë. 2020 ã.
 *      Author: dim
 */

#ifndef SRC_ARC_DMA_DRIVER_H_
#define SRC_ARC_DMA_DRIVER_H_


#define DRIVER_IOCTL_MAGIC              'r'
#define DRIVER_IOCTL_START_STREAM       _IOW ( DRIVER_IOCTL_MAGIC, 3, int )
#define DRIVER_IOCTL_READ_STATE         _IOR ( DRIVER_IOCTL_MAGIC, 4, int )
#define DRIVER_IOCTL_WRITE_STATE        _IOW ( DRIVER_IOCTL_MAGIC, 5, int )


#define FFT_CORE_DMA_BUFFER_SIZE          0x01400000 // 20 Mb

int fd_fft_0;
int fd_fft_1;
unsigned int *fft_0_core_buffer;
unsigned int *fft_1_core_buffer;

unsigned int actual_bytes;

int  fft_0_dma_init(int N_packets);
void fft_0_dma_start(int N_packets);
int  fft_0_dma_stop();
void fft_0_dma_exit();


int  fft_1_dma_init(int N_packets);
void fft_1_dma_start(int N_packets);
int  fft_1_dma_stop();
void fft_1_dma_exit();



#endif /* SRC_ARC_DMA_DRIVER_H_ */
