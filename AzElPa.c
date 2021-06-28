#include "AzElPa.h"
#include <math.h>
#include "iodkf.h"
/**
* @file AzElPa.c
*
* @brief  Implementación de la funcion: void AzElPa(double* v, double *Az, double *El, double* dAds, double* dEds)
*
* @author Juan Pérez Garrido
*
*/
void AzElPa(double* v, double *Az, double *El, double* dAds, double* dEds) {
	//como no podemos devolver todos los parámetros, los ponemos como entrada/salida
	double rho = sqrt(v[0] * v[0] + v[1] * v[1]);

	*Az = atan2(v[0], v[1]);
	if (*Az < 0.0) {
		*Az = *Az + pi2;
	}
	*El = atan(v[2] / rho);

	dAds[0] = v[1] / (rho * rho);
	dAds[1] = -v[0] / (rho * rho);
	dAds[2] = 0.0;

	double d = dot(v, 3, v, 3);
	
	dEds[0] = (-v[0] * v[2] / rho);
	dEds[1] = (-v[1] * v[2] / rho);
	dEds[2] = rho;
	double* temp = vectorPorEscalar(dEds, 3, (1.0 / d));
	dEds[0] = temp[0];
	dEds[1] = temp[1];
	dEds[2] = temp[2];

}