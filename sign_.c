#pragma once
#include "sign_.h"
/**
* @file sign_.c
*
* @brief  Implementación de la funcion: double sign_(double a, double b) 
*
* @author Juan Pérez Garrido
*
*/
double sign_(double a, double b) {
	double resultado;
	if (b >= 0.0) {
		resultado = fabs(a);
	}
	else {
		resultado = -fabs(a);
	}
	return resultado;
}