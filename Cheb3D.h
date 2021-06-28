#pragma once
double* Cheb3D(double t,double N,double Ta,double Tb,
	double *Cx,double *Cy,double *Cz);

/**
* @file Cheb3D.h
* @brief Es la aproximación de Chebysev de vectores tridimensionales
*
* @param t
* @param N Número de coeficientes
* @param Ta Inicio del intervalo
* @param Tb Final del intervalo
* @param Cx Coeficientes del polinomio de Chebyshev (coordenada x)
* @param Cy Coeficientes del polinomio de Chebyshev (coordenada y)
* @param Cz Coeficientes del polinomio de Chebyshev (coordenada z)
* @return Devuelve el vector con la aproximación de Chebyshev de vectores tridimensionales.
*
*/