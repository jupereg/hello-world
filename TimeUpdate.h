#pragma once
double** TimeUpdate(double** P,int f1,int c1, double** Phi,int f2,
	int c2, double Qdt);
/**
* @file TimeUpdate.h
* @param P Matriz
* @param f1 n�mero de filas de P
* @param c1 n�mero de columnas de P
* @param Phi Matriz
* @param f2 n�mero de filas de Phi
* @param c2 n�mero de columnas de Phi
* @param Qdt
*
* @return Devuelve otra matriz
*
*/