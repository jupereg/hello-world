#include "Cheb3D.h"
#include<stdlib.h>
#include<math.h>
#include <stdio.h>
#include "vector.h"
#include "matriz.h"
/**
* @file Cheb3D.c
*
* @brief  Implementación de la función:
*	double* Cheb3D(double t, double N, double Ta, double Tb, double* Cx, double* Cy, double* Cz)
*
* @author Juan Pérez Garrido
*
*/

double* Cheb3D(double t, double N, double Ta, double Tb,
	double* Cx, double* Cy, double* Cz) {

	double** old_f1 = NULL, ** t1 = NULL, 
		** t2 = NULL, ** t3 = NULL, ** t4 = NULL;


	if ((t < Ta) || (Tb < t))
	{
		printf("ERROR: Time out of range in Cheb3D::Value\n");
		exit(EXIT_FAILURE);
	}
	//Clenshaw algorithm
	double tau = (2 * t - Ta - Tb) / (Tb - Ta);

	double **f1 = crearMatriz(1, 3);
	double **f2 = crearMatriz(1, 3);

	for (int i = N; i >= 2; i--)
	{
		old_f1 = f1;

		t1 = matrizPorEscalar(f1, 1, 3, 2 * tau);
		t2 = matrizPorEscalar(f2, 1, 3, -1.0);
		t3 = sumaMatrices(t1, 1, 3, t2, 1, 3);

		t4 = crearMatriz(1, 3);
		t4[0][0] = Cx[i - 1]; t4[0][1] = Cy[i - 1]; t4[0][2] = Cz[i - 1];

		f1 = sumaMatrices(t3, 1, 3, t4, 1, 3);

		f2 = old_f1;
	}

	double **t5 = matrizPorEscalar(f1, 1, 3, tau);
	double **t6 = matrizPorEscalar(f2, 1, 3, -1.0);
	double **t7 = sumaMatrices(t5, 1, 3, t6, 1, 3);

	double **t8 = crearMatriz(1, 3);
	t8[0][0] = Cx[0];
	t8[0][1] = Cy[0];
	t8[0][2] = Cz[0];

	double** ChebApp = sumaMatrices(t7, 1, 3, t8, 1, 3);

	liberarMatriz(old_f1, 1, 3);
	liberarMatriz(f1, 1, 3);
	liberarMatriz(t1, 1, 3);
	liberarMatriz(t2, 1, 3);
	liberarMatriz(t3, 1, 3);
	liberarMatriz(t4, 1, 3);
	liberarMatriz(t5, 1, 3);
	liberarMatriz(t6, 1, 3);
	liberarMatriz(t7, 1, 3);
	liberarMatriz(t8, 1, 3);


	return ChebApp;
}