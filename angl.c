#include "angl.h"
#include <math.h>
/**
* @file angl.c
*
* @brief  Implementación de la funcion: double angl(double* v1,int tam1, double* v2,int tam2)
*
* @author Juan Pérez Garrido
*
*/
double angl(double* v1,int tam1, double* v2,int tam2) {
	double small = 0.00000001;
	double undefinded = 999999.1;
	double theta;

	double magv1 = normaVector(v1, tam1);
	double magv2 = normaVector(v2, tam2);

	if (magv1 * magv2 > small^2) {
		double temp = dot(v1, tam1, v2, tam2)/(magv1*magv2);
		if (fabs(temp) > 1.0) {
			if (temp > 0)
				temp = 1 * 1.0;
			if (temp == 0)
				temp = 0 * 1.0;
			if (temp < 0)
				temp = -1 * 1.0;
		}
		theta = acos(temp);
	}
	else {
		theta = undefinded;
	}
	return theta;
}