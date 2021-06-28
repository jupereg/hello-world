#include "matriz.h"
/**
* @file vector.c
*
* @brief  Implementación de todas las operaciones con vectores.
*
* @author Juan Pérez Garrido
*
*/
double* crearVector(int tam) {
	double* vector;
	vector = (double*)malloc(tam * sizeof(double));
	for (int i = 0; i < tam; i++) {
		vector[i] = 1.0;
	}
	return vector;
}
void imprimirVector(double* v, int tam) {
	for (int i = 0; i < tam; i++) {
		printf("%f \n", v[i]);
	}
}
void liberarVector(double* v) {
	free(v);
}
int compararVectores(double* v1, double* v2, int tam1, int tam2) {
	int i, j;
	if (tam1 != tam2) {
		printf("vec: different dimension\n");
		exit(EXIT_FAILURE);
		return 0;
	}
	for (int i = 0; i < tam1; i++) {
		if (fabs(v1[i] - v2[i]) > pow(10, -6))
			return 0;
		return 1;
	}
}
double* vectorPorEscalar(double* v, int tam, double escalar) {
	double* v1 = crearVector(tam);
	for (int i = 0; i < tam; i++) {
		v1[i] = v[i] * escalar;
	}
	return v1;
}
double normaVector(double* v, int tam) {
	double suma = 0.0;
	double resultado;
	for (int i = 0; i < tam; i++) {
		suma = suma + v[i] * v[i];
	}
	resultado = sqrt(suma);
}
double* vectorMasEscalar(double* v, int tam, double escalar) {
	double* sol = crearVector(tam);
	for (int i = 0; i < tam; i++) {
		sol[i] = v[i] + escalar;
	}
	return sol;
}

//para JPL_EPh
int find1(double* v1, int n1, double* v2, int n2, double JD) {
	int i = 0,ind=0;
	if (n1 != n2) {
		printf("vec: different dimension\n");
		exit(EXIT_FAILURE);
		return 0;
	}
	while (i < n1) {
		if (v1[i] <= JD) {
			ind = i;
		}
		i++;
	}
	if (ind == 0) {
		i = 0;
		while (i < n2) {
			if (v2[i] > JD) {
				ind = i;
			}
			i++;
		}
	}
	return ind;
}
double dot(double* v1, int tam1, double* v2, int tam2) {
	if (tam1 != tam2) {
		exit(EXIT_FAILURE);
	}
	double resultado = 0;
	for (int i = 0; i < tam1; i++) {
		resultado = resultado + (v1[i] * v2[i]);
	}
	return resultado;
}

double* dotVectorial(double* v1, int tam1, double* v2, int tam2) {
	//en este caso devolveremos el vector resultado del producto vectorial
	//solo para vectores de 3 componentes
	if (tam1 != tam2 || tam1 != 3 || tam2 != 3) {
		exit(EXIT_FAILURE);
	}
	double* resultado = crearVector(tam1);
	resultado[0] = v1[1] * v2[2] - v1[2] * v2[1];
	resultado[1] = v1[2] * v2[0] - v1[0] * v2[2];
	resultado[2] = v1[0] * v2[1] - v1[1] * v2[0];
	return resultado;
}
double* sumaVectores(double* v1, int t1, double* v2, int t2) {
	//el tamaño tiene que ser igual
	double* resultado = crearVector(t1);
	if (t1 == t2) {
		for (int i = 0; i < t1; i++) {
			resultado[i] = v1[i] + v2[i];
		}
	}
	return resultado;
}
double* restaVectores(double* v1, int t1, double* v2, int t2) {
	//el tamaño tiene que ser igual
	double* resultado = crearVector(t1);
	if (t1 == t2) {
		for (int i = 0; i < t1; i++) {
			resultado[i] = v1[i] - v2[i];
		}
	}
	return resultado;
}
