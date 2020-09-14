#ifndef SRC_RECOGNITION_H_
#define SRC_RECOGNITION_H_

typedef enum {
	Unclassified,
	Non_radar_signal,
	Iskra,
	MultaNova,
	Multa_CD,
	MultaNova6F,
	Ramera_222,
	Stalker_34G
} en_Types;


extern char *recognition (st_TracksInfo *TracksInfo, uint *TracksLengths, uint Ntraces, uint Nstep, st_trace **traces, uint Nwindow, float *Prob_of_recogn);

#endif /* SRC_RECOGNITION_H_ */
