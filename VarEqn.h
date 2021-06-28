#pragma once
#include "IERS.h"
#include "matriz.h"
#include "vector.h"
#include "PrecMatrix.h"
#include "NutMat.h"
#include "PoleMatrix.h"
#include "iodkf.h"
#include "GHAMatrix.h"
#include "AccelHarmonic.h"
#include "G_AccelHarmonic.h"
#include "timediff.h"
#include "globales.h"
void VarEqn(double x, double* yPhi, double* yPhip);
/**
* @file VarEqn
* @brief Calcula las ecuaciones variacionales, es decir, la derivada del  vector de estado y la matriz de transición de estado.
* 
*
* @param x Tiempo desde época en [s]
* @param yPhi (6 + 36) -dim vector que comprende el vector de estado(y) y la matriz de transición de estado (Phi) en orden de almacenamiento por columnas.
* @param yPhip Devuelve la derivada de yPhi
*
*/
