#pragma once
double* AccelPointMass(double* r, int nr, double* s, int ns,
	double GM);

/**
*
* @file AccelPointMass.h
* @brief Calcula la aceleración por perturbación debida a una masa puntual.
*
* @param r Posición del vector satelital
* @param nr Tamaño de r
* @param s Point mass position vector
* @param ns tamaño de s
* @param GM Coeficiente gravitacional del punto de masa
* @return Vector de aceleración --> (a=d^2r/dt^2)
*
*/