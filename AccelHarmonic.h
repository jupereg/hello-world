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
* @brief Calcula la aceleración debida al campo gravitatorio armónico del cuerpo central
*
* @param r Vector de posición del satélite en el sistema inercial
* @param E Matriz de transformación a sistema de cuerpo fijo
* @param n_max Grado máximo
* @param m_max Orden máxima (m_max<=n_max; m_max=0 for zonals, only)
*
* @return Vector de aceleración --> (a=d^2r/dt^2)
*
*/