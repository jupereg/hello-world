#include "TimeUpdate.h"
#include <stdlib.h>
#include "matriz.h"
/**
* @file TimeUpdate.c
*
* @brief  Implementación de la funcion: double** TimeUpdate(double** P, int f1, int c1, double** Phi, int f2, int c2, double Qdt)
*
* @author Juan Pérez Garrido
*
*/
double** TimeUpdate(double** P, int f1, int c1, double** Phi, int f2, int c2, double Qdt) {
	if (f1 != f2 || c1 != c2)
		exit(EXIT_FAILURE);

	P = productoMatrices(Phi, f2, c2, P, f1, c1);
	double** Phit = crearMatriz(f2, c2);
	Phit = matrizTraspuesta(Phi, f2);
	P = productoMatrices(P, f1, c1, Phit, f2, c2);
	for (int i = 0; i < f1; i++) {
		for (int j = 0; j < f1; j++) {
			P[i][j] = P[i][j] + Qdt;
		}
	}
	return P;
}