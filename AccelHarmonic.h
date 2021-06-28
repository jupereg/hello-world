#pragma once
#include <math.h>
#include "Legendre.h"
#include "iodkf.h"
#include "matriz.h"
#include "vector.h"
#include "globales.h"
double* AccelHarmonic(double* r, double** E, int n_max, int m_max);


/**
*
* @file AccelHarmonic.h
* @brief Calcula la aceleraci�n debida al campo gravitatorio arm�nico del cuerpo central
*
* @param r Vector de posici�n del sat�lite en el sistema inercial
* @param E Matriz de transformaci�n a sistema de cuerpo fijo
* @param n_max Grado m�ximo
* @param m_max Orden m�xima (m_max<=n_max; m_max=0 for zonals, only)
*
* @return Vector de aceleraci�n --> (a=d^2r/dt^2)
*
*/