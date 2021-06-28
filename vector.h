#pragma once
#include <stdio.h>
double* crearVector(int tam);
void imprimirVector(double* v, int tam);
void liberarVector(double* v);
int compararVectores(double* v1, double* v2, int tam1, int tam2);
double* vectorPorEscalar(double* v, int tam, double escalar);
double* vectorMasEscalar(double* v, int tam, double escalar);
double normaVector(double* v, int tam);
int find1(double* v1, int n1, double* v2, int n2, double JD);
double dot(double* v1, int tam1, double* v2, int tam2);
double* dotVectorial(double* v1, int tam1, double* v2, int tam2);
double* sumaVectores(double* v1, int t1, double* v2, int t2);
double* restaVectores(double* v1, int t1, double* v2, int t2);

/**
* @file vector.h
* @brief Incluye las operaciones con vectores
*/