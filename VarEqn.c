#include "VarEqn.h"

/**
* @file VarEqn.c
*
* @brief  Contiene la implementacion de la funcion void VarEqn(double x, double* yPhi, double *yPhip)
*
* @author Juan Pérez Garrido
*
*/

void VarEqn(double x, double* yPhi, double* yPhip)
{
	extern Param AuxParam;
	extern double** eopdata,**Cnm,**Snm;

	double x_pole = 0, y_pole = 0, UT1_UTC = 0, LOD = 0, dpsi = 0, deps = 0, dx_pole = 0, dy_pole = 0, TAI_UTC = 0, Mjd_UT1;
	double UT1_TAI = 0, UTC_GPS = 0, UT1_GPS = 0, TT_UTC = 0, GPS_UTC = 0, * v, * r, ** Phi;

	IERS(AuxParam.Mjd_UTC, "l", &x_pole, &y_pole, &UT1_UTC, &LOD, &dpsi, &deps, &dx_pole, &dy_pole, &TAI_UTC);

	timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);

	Mjd_UT1 = AuxParam.Mjd_TT + (UT1_UTC - TT_UTC) / 86400;

	/* Matriz de transformacion */

	double** P = PrecMat(MJD_J2000, AuxParam.Mjd_TT + x / 86400);
	double** N = NutMat(AuxParam.Mjd_TT + x / 86400);

	double** T = productoMatrices(N, 3, 3, P, 3, 3);
	double** E = productoMatrices(productoMatrices(PoleMatrix(x_pole, y_pole), 3, 3, GHAMatrix(Mjd_UT1), 3, 3), 3, 3, T, 3, 3);

	/* Componentes del vector de estado */

	r = crearVector(3);
	r[0] = yPhi[0]; r[1] = yPhi[1]; r[2] = yPhi[2];

	v = crearVector(3);
	v[0] = yPhi[3]; v[1] = yPhi[4]; v[2] = yPhi[5];

	Phi = crearMatriz(6, 6);

	/* Matriz de transición de estado */
	for (int j = 1; j <= 6; j++)
	{
		//Phi(:, j) = yPhi(6 * j + 1:6 * j + 6);

		int incio = 6 * j + 1;
		int fin = 6 * j + 6;
		int acum = incio;

		double* aux = crearVector(fin - incio + 1);
		for (int i = 0; i < (fin - incio + 1); i++)
		{
			aux[i] = yPhi[acum - 1];
			acum++;
		}

		for (int i = 0; i < 6; i++)
		{
			Phi[i][j - 1] = aux[i];
		}

		liberarVector(aux);

	}

	/* Aceleración y gradiente */

	double* a = AccelHarmonic(r, E, AuxParam.n, AuxParam.m);
	double** G = G_AccelHarmonic(r, E, AuxParam.n, AuxParam.m);

	/* Derivado temporal de la matriz de transición de estado */

	double** dfdy = crearMatriz(6, 6);

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			dfdy[i][j] = 0.0;
			dfdy[i + 3][j] = G[i][j];
			if (i == j)
			{
				dfdy[i][j + 3] = 1;
			}
			else
			{
				dfdy[i][j + 3] = 0;
			}
			dfdy[i + 3][j + 3] = 0.0;
		}
	}

	double** Phip = productoMatrices(dfdy, 6, 6, Phi, 6, 6);

	/* Derivada del vector de estado combinado y la matriz de transición de estado */

	for (int i = 0; i < 3; i++)
	{
		yPhip[i] = v[i];
		yPhip[i + 3] = a[i];
	}

	for (int i = 1; i <= 6; i++)
	{
		for (int j = 1; j <= 6; j++)
		{
			yPhip[6 * j + i - 1] = Phip[i - 1][j - 1];
		}
	}

	liberarVector(v);
	liberarVector(r);
	liberarMatriz(Phi, 6, 6);

	liberarMatriz(P, 3, 3);
	liberarMatriz(N, 3, 3);
	liberarMatriz(T, 3, 3);
	liberarMatriz(E, 3, 3);

	liberarVector(a);
	liberarMatriz(G, 3, 3);
	liberarMatriz(dfdy, 6, 6);
	liberarMatriz(Phip, 6, 6);


}