#pragma once
double* AccelPointMass(double* r, int nr, double* s, int ns,
	double GM);

/**
*
* @file AccelPointMass.h
* @brief Calcula la aceleraci�n por perturbaci�n debida a una masa puntual.
*
* @param r Posici�n del vector satelital
* @param nr Tama�o de r
* @param s Point mass position vector
* @param ns tama�o de s
* @param GM Coeficiente gravitacional del punto de masa
* @return Vector de aceleraci�n --> (a=d^2r/dt^2)
*
*/