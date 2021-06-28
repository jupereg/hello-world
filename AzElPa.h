#pragma once
#include "vector.h"
void AzElPa(double* v,double *Az,double *El,double *dAds,double *dEds);
/**
* @file AzElPa.h
* @brief Calcula azimuth, la elevaci�n y los parciales a partir de las coordenadas de la tangente local.
*
* @param s vector de coordenadas de la tangente local topoc�ntrica 
* @param Az Azimuth [rad]
* @param El Elevacion [rad]
* @param dAds Parciales of azimuth w.r.t. s
* @param dEds Parciales de elevacion w.r.t. s
*/