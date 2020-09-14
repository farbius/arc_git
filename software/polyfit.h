#ifndef _POLIFITGSL_H
#define _POLIFITGSL_H
#include "gsl/gsl_multifit.h"
#include <stdbool.h>
#include <math.h>
bool polyfit(int obs, int degree,
		   double *dx, double *dy, double *store); /* n, p */
#endif



/*#ifndef SRC_POLYFIT_H_
#define SRC_POLYFIT_H_


extern int polyfit(const double* const dependentValues,
        const double* const independentValues,
        unsigned int        countOfElements,
        unsigned int        order,
        double*             coefficients);


//#endif /* SRC_POLYFIT_H_ */
