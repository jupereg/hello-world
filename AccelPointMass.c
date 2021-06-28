#include "AccelPointMass.h"
#include "vector.h"
#include <math.h>
/**
* @file AccelPointMass.c
*
* @brief  Implementaci�n de la funcion: double* AccelPointMass(double* r, int nr, double* s, int ns, double GM)
*
* @author Juan P�rez Garrido
*
*/
double* AccelPointMass(double* r, int nr, double* s, int ns, double GM) {
	
	double* sum = sumaVectores(r, nr, vectorPorEscalar(s, ns, -1.0), ns);

	return sumaVectores(vectorPorEscalar(sum,ns,-GM/pow(normaVector(sum,3),3.0)),ns, 
		vectorPorEscalar(s, ns, -GM / pow(normaVector(s, 3), 3.0),ns),ns);
		
}