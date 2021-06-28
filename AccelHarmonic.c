#include "AccelHarmonic.h"
/**
* @file Legendre.c
*
* @brief  Implementación de la funcion: AccelHarmonic(double* r, double** E, int n_max, int m_max)
*
* @author Juan Pérez Garrido
*
*/
double* AccelHarmonic(double* r, double** E, int n_max, int m_max)
{
	double  b1, b2, b3;
	extern double** Cnm,**Snm;

	double r_ref = 6378.1363e3;
	double gm = 398600.4415e9;


	double **aux_r_bf = matrizPorVector(E, 3, 3, r, 3);
	double *r_bf = crearVector(3);
	r_bf[0] = aux_r_bf[0][0]; r_bf[1] = aux_r_bf[1][0]; r_bf[2] = aux_r_bf[2][0];

	// Auxiliary quantities
	double d = normaVector(r_bf, 3); 
	double latgc = asin(r_bf[2] / d);
	double lon = atan2(r_bf[1], r_bf[0]);

	double **pnm = crearMatriz(n_max + 1, m_max + 1);
	double **dpnm = crearMatriz(n_max + 1, m_max + 1);

	Legendre(n_max, m_max, latgc, pnm, dpnm);

	double dUdr = 0;
	double dUdlatgc = 0;
	double dUdlon = 0;
	double q3 = 0; 
	double q2 = q3; 
	double q1 = q2;

	for (int n = 0; n <= n_max; n++)
	{
		b1 = ((-1.0 * gm) / (d * d)) * pow((r_ref / d), n) * (n + 1.0);
		b2 = (gm / d) * pow((r_ref / d), n);
		b3 = (gm / d) * pow((r_ref / d), n);

		for (int m = 0; m <= m_max; m++)
		{
			q1 = q1 + pnm[n][m] * (Cnm[n][m] * cos(m * lon) + Snm[n][m] * sin(m * lon));
			q2 = q2 + dpnm[n][m] * (Cnm[n][m] * cos(m * lon) + Snm[n][m] * sin(m * lon));
			q3 = q3 + m * pnm[n][m] * (Snm[n][m] * cos(m * lon) - Cnm[n][m] * sin(m * lon));
		}

		dUdr = dUdr + q1 * b1;
		dUdlatgc = dUdlatgc + q2 * b2;
		dUdlon = dUdlon + q3 * b3;
		q3 = 0; q2 = q3; q1 = q2;
	}

	// Body-fixed position 

	double r2xy = r_bf[0] * r_bf[0] + r_bf[1] * r_bf[1];

	double ax = (1 / d * dUdr - r_bf[2] / (d * d * sqrt(r2xy)) * dUdlatgc) * r_bf[0] - (1 / r2xy * dUdlon) * r_bf[1];
	double ay = (1 / d * dUdr - r_bf[2] / (d * d * sqrt(r2xy)) * dUdlatgc) * r_bf[1] + (1 / r2xy * dUdlon) * r_bf[0];
	double az = 1 / d * dUdr * r_bf[2] + sqrt(r2xy) / (d * d) * dUdlatgc;

	double *a_bf = crearVector(3);
	a_bf[0] = ax; a_bf[1] = ay; a_bf[2] = az;

	//Inertial acceleration

	double **aux = matrizPorVector(matrizTraspuesta(E, 3, 3), 3, 3, a_bf, 3);
	double *a = crearVector(3);
	a[0] = aux[0][0]; a[1] = aux[1][0]; a[2] = aux[2][0];

	liberarMatriz(aux_r_bf, 3, 3);
	liberarVector(r_bf);
	liberarMatriz(pnm, n_max + 1, m_max + 1);
	liberarMatriz(dpnm, n_max + 1, m_max + 1);
	liberarVector(a_bf);
	liberarMatriz(aux, 1, 3);

	return a;
}