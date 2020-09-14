#include "stdio.h"
#include <stdlib.h>
#include "Soft_Part_I.h"


void Soft_Part_I (st_Results *R, float *Sf, uint *Begs, uint *Fins, uint Width) {

	// 0. variables:
	uint 	numoflines;
	uchar 	needdeleting;
	uint 	Nlength 	= 1024;
	uint 	n_iter 		= 0;
	uint 	*Lens 		= malloc(sizeof(uint)*1024);
	float 	*Bets 		= malloc(sizeof(float)*1024);
	uchar 	*IsJoin 	= malloc(sizeof(uchar)*1024);
	float 	Spare 		= 5.0;
	uint 	i, j, x;
	float 	Center;
	float 	P;

	// 1. check Begs:
	if (Begs[0] == 0) {
		R->Begs[0] 		= 0;
		R->Fins[0] 		= 0;
		R->Types[0] 	= 0;
		R->Centers[0] 	= 0;
		R->Widths[0] 	= 0;
		R->Num 			= 0;
		R->n_iter 		= 0;
		return;
	}

	// 2. processing:
	while (1) {

		// 2.1. lengths of lines:
		for (numoflines = 0; numoflines < 1024; numoflines++) {
			if (Begs[numoflines] == 0) {
				Lens[numoflines] = 0;
				break;
			} else {
				Lens[numoflines] = Fins[numoflines] - Begs[numoflines];
			}
		}

		// 2.2. check number of lines:
		if (numoflines == 1) break;

		// 2.3. distances between lines & flag of lines concatenation:
		for (i = 0; i < numoflines; i++) {
			Bets[i] 	= 0;
			IsJoin[i] 	= 0;
		}

		// 2.4. check for a concatenation:
		needdeleting = 0;
		for (i = 0; i < numoflines-1; i++) {
			Bets[i] = Spare * (Begs[i+1] - Fins[i]);
			if ((Begs[i+1] > 0) && (Fins[i] > 0)) {
				if ((Lens[i] >= Bets[i]) || (Bets[i] <= Lens[i+1])) {
					IsJoin[i] = 1; 										// mark of space deleting
					needdeleting = 1;
				}
			}
		}

		if (needdeleting == 0) break;

		// 2.5. concatination:
		j = 0;
		for (i = 0; i < numoflines; i++) {
			if (IsJoin[i] == 0) {
				Begs[j+1] 	= Begs[i+1];
				Fins[j] 	= Fins[i];
				j++;
			}
		}
		numoflines = j;
		Begs[numoflines] = 0;
		Fins[numoflines] = 0;

		// 2.6. check of a iterations:
		n_iter++;
		if (n_iter > 50) {
			fprintf (stderr, "Soft_Part_I: The number of iterations exceeds 50!!!\n");
			break;
		}

	}

	// 3. Result:
	for (i = 0; i < numoflines; i++) {

		// 3.1. copy lines:
		R->Begs[i] = Begs[i];
		R->Fins[i] = Fins[i];

		// 3.2. calc of center and width:
		Center 	= 0;
		P 		= 0;
		for (x = Begs[i]; x <= Fins[i]; x++) {
			Center 	+= Sf[x-1] * Sf[x-1] * x;
			P 		+= Sf[x-1] * Sf[x-1];
		}
		R->Centers[i] 	= Center / P;
		R->Widths[i] 	= R->Fins[i] - R->Begs[i] + 1;
		if (R->Widths[i] > Nlength) R->Widths[i] = Nlength; 				// ??? I don't understand that this

		// 3.3. other data:
		R->Types[i] = (R->Widths[i] > Width);

	}
	R->Num 		= numoflines;
	R->n_iter 	= n_iter;
	R->Begs[i] 	= 0;
	R->Fins[i] 	= 0;

	// 4. cleaning:
	free(Lens);
	free(Bets);
	free(IsJoin);

}
