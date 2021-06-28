#include "Accel.h"
#include "globales.h"

/**
* @file Accel.c
*
* @brief  Implementación de la funcion: void Accel(double x, double* Y, double* dY)
*
* @author Juan Pérez Garrido
*
*/

void Accel(double x, double* Y, double* dY)
{
	extern Param AuxParam;
	extern double** eopdata,**PC,**Cnm,**Snm;


	double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
	IERS(AuxParam.Mjd_UTC + x / 86400, "l", &x_pole, &y_pole, &UT1_UTC, &LOD, &dpsi, &deps, &dx_pole, &dy_pole, &TAI_UTC);

	double  UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
	timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);

	double Mjd_UT1 = AuxParam.Mjd_UTC + x / 86400 + UT1_UTC / 86400;
	double Mjd_TT = AuxParam.Mjd_UTC + x / 86400 + TT_UTC / 86400;

	double** P = PrecMat(MJD_J2000, Mjd_TT);
	double** N = NutMat(Mjd_TT);
	double** T = productoMatrices(N, 3, 3, P, 3, 3);
	double** E = productoMatrices(productoMatrices(PoleMatrix(x_pole, y_pole), 3, 3, GHAMatrix(Mjd_UT1), 3, 3), 3, 3, T, 3, 3);


	double MJD_TDB = Mjday_TDB(Mjd_TT);

	double* r_Mercury = crearVector(3), * r_Venus = crearVector(3), * r_Earth = crearVector(3), * r_Mars = crearVector(3), * r_Jupiter = crearVector(3), * r_Saturn = crearVector(3), * r_Uranus = crearVector(3), * r_Neptune = crearVector(3), * r_Pluto = crearVector(3), * r_Moon = crearVector(3), * r_Sun = crearVector(3);
	JPL_Eph_DE430(MJD_TDB, r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun);


	//Acceleration due to harmonic gravity field
	double* aux = crearVector(3);
	aux[0] = Y[0]; aux[1] = Y[1]; aux[2] = Y[2];

	double* a = AccelHarmonic(aux, E, AuxParam.n, AuxParam.m);

	//Luni-solar perturbations
	if (AuxParam.sun)
	{
		a = sumaVectores(a, 3, AccelPointMass(aux, 3, r_Sun, 3, GM_Sun), 3);
	}

	if (AuxParam.moon)
	{
		a = sumaVectores(a, 3, AccelPointMass(aux, 3, r_Moon, 3, GM_Moon), 3);
	}

	//Planetary perturbations
	if (AuxParam.planets)
	{
		a = sumaVectores(a, 3, AccelPointMass(aux, 3, r_Mercury, 3, GM_Mercury), 3);
		a = sumaVectores(a, 3, AccelPointMass(aux, 3, r_Venus, 3, GM_Venus), 3);
		a = sumaVectores(a, 3, AccelPointMass(aux, 3, r_Mars, 3, GM_Mars), 3);
		a = sumaVectores(a, 3, AccelPointMass(aux, 3, r_Jupiter, 3, GM_Jupiter), 3);
		a = sumaVectores(a, 3, AccelPointMass(aux, 3, r_Saturn, 3, GM_Saturn), 3);
		a = sumaVectores(a, 3, AccelPointMass(aux, 3, r_Uranus, 3, GM_Uranus), 3);
		a = sumaVectores(a, 3, AccelPointMass(aux, 3, r_Neptune, 3, GM_Neptune), 3);
		a = sumaVectores(a, 3, AccelPointMass(aux, 3, r_Pluto, 3, GM_Pluto), 3);
	}


	dY[0] = Y[3]; 
	dY[1] = Y[4]; 
	dY[2] = Y[5];
	dY[3] = a[0]; 
	dY[4] = a[1]; 
	dY[5] = a[2];

	liberarMatriz(P, 3, 3);liberarMatriz(N, 3, 3);
	liberarMatriz(T, 3, 3);liberarMatriz(E, 3, 3);

	liberarVector(r_Mercury);liberarVector(r_Venus);liberarVector(r_Earth);
	liberarVector(r_Mars);liberarVector(r_Jupiter);liberarVector(r_Saturn);
	liberarVector(r_Uranus);liberarVector(r_Neptune);liberarVector(r_Pluto);
	liberarVector(r_Moon);liberarVector(r_Sun);liberarVector(aux);liberarVector(a);

}