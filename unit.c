#include "unit.h"
/**
* @file unit.c
*
* @brief  Implementación de la funcion: double* Unit(double* v,int tam)
*
* @author Juan Pérez Garrido
*
*/
double* Unit(double* v,int tam) {
	double* outvec = crearVector(3);
	double small = 0.000001;
	double magv = normaVector(v,tam);

	if (magv > small) {
		for (int i = 0; i < 3; i++)
			outvec[i] = v[i] / magv;
	}
	else {
		for (int i = 0; i < 3; i++) {
			outvec[i] = 0.0;
		}
	}
	return outvec;
}