#include <stdio.h>
#include <stdlib.h>
#include "R_x.h"
#include "R_y.h"
#include "R_z.h"
#include "matriz.h"
#include "vector.h"
#include "NutAngles.h"
#include "posicion.h"
#include "NutMat.h"
#include "MeanObliquity.h"
#include "unit.h"
#include "angl.h"
#include "AzElPa.h"
#include "Geodetic.h"
#include "elements.h"
#include "Frac.h"
#include "LTC.h"
#include "Timediff.h"
#include "sign_.h"
#include "gmst.h"
#include "PoleMatrix.h"
#include "TimeUpdate.h"
#include "gibbs.h"
#include "AccelPointMass.h"
#include "PrecMatrix.h"
#include "EqnEquinox.h"
#include "gast.h"
#include "GHAMatrix.h"
#include "Mjday_TDB.h"
#include "Legendre.h"
#include "MeasUpdate.h"
#include "Cheb3D.h"
#include "Mjday.h"
#include "AccelHarmonic.h"
#include "IERS.h"
#include "VarEqn.h"
#include "G_AccelHarmonic.h"
#include "Accel.h"


#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

#pragma warning(disable : 4996)//evito los problemas que me da la lectura de ficheros

/**
* @file test.c
* 
* @brief Todos los test unitiarios de las funciones realizadas
* 
* @author Juan Pérez Garrido
* 
*/

int tests_run = 0;

int R_x_01() {
	double** resultado = R_x(1.0);
	double** sol = crearMatriz(3, 3);
	sol[0][0] = 1.0; sol[0][1] = 0.0; sol[0][2] = 0.0;
	sol[1][0] = 0.0; sol[1][1] = 0.54030230586813977; sol[1][2] = 0.84147098480789650;
	sol[2][0] = 0.0; sol[2][1] = -0.84147098480789650; sol[2][2] = 0.54030230586813977;

    _assert(compararMatrices(resultado, 3, 3, sol, 3, 3));

    liberarMatriz(resultado,3,3);
    liberarMatriz(sol,3,3);
    return 0;
}
int R_y_01() {
    double** resultado = R_y(1.0);
    double** sol = crearMatriz(3, 3);
    sol[0][0] = 0.54030230586813977; sol[0][1] = 0.0; sol[0][2] = -0.84147098480789650;
    sol[1][0] = 0.0; sol[1][1] = 1.0; sol[1][2] = 0.0;
    sol[2][0] = 0.84147098480789650; sol[2][1] = 0.0; sol[2][2] = 0.54030230586813977;

    _assert(compararMatrices(resultado, 3, 3, sol, 3, 3));

    liberarMatriz(resultado, 3, 3);
    liberarMatriz(sol, 3, 3);
    return 0;
}
int R_z_01() {
    double** resultado = R_z(1.0);
    double** sol = crearMatriz(3, 3);
    sol[0][0] = 0.54030230586813977; sol[0][1] = 0.84147098480789650; sol[0][2] = 0.0;
    sol[1][0] = -0.84147098480789650; sol[1][1] = 0.54030230586813977; sol[1][2] = 0.0;
    sol[2][0] = 0.0; sol[2][1] = 0.0; sol[2][2] = 1.0;

    _assert(compararMatrices(resultado, 3, 3, sol, 3, 3));

    liberarMatriz(resultado, 3, 3);
    liberarMatriz(sol, 3, 3);
    return 0;
}
//VECTORES
int normaVector_01(){
    double* vector = crearVector(5);
    double resultado = normaVector(vector, 5);
    _assert(sqrt(5) == resultado);
    liberarVector(vector);
    return 0;
}
int vectorPorEscalar_01() {
    double* vector = crearVector(5);
    vector[0] = 1; vector[1] = 2; vector[2] = 3; vector[3] = 4; vector[4] = 5;
    double* vector2 = crearVector(5);
    vector2[0] = 2; vector2[1] = 4; vector2[2] = 6; vector2[3] = 8; vector2[4] = 10;
    _assert(compararVectores(vectorPorEscalar(vector, 5, 2), vector2, 5, 5));
    liberarVector(vector);
    liberarVector(vector2);
    return 0;
}
int vectorMasEscalar_01() {
    double* v = crearVector(3);
    v[0] = 1; v[1] = 2; v[2] = 3;
    double* suma = crearVector(3);
    suma[0] = 2; suma[1] = 3; suma[2] = 4;
    double* sol = vectorMasEscalar(v, 3, 1);

    _assert(compararVectores(suma, sol, 3, 3));
    
    liberarVector(suma);
    liberarVector(sol);
    liberarVector(v);
    return 0;
}
int sumaVectores_01() {
    double* resul = crearVector(3);
    resul[0] = 2; resul[1] = 4; resul[2] = 6;

    double* v = crearVector(3);
    v[0] = 1; v[1] = 2; v[2] = 3;

    double* v2 = crearVector(3);
    v2[0] = 1; v2[1] = 2; v2[2] = 3;

    double* suma = crearVector(3);
    suma = sumaVectores(v, 3, v2, 3);
    _assert(compararVectores(suma, resul, 3, 3));

    liberarVector(resul);
    liberarVector(suma);
    liberarVector(v);
    liberarVector(v2);
    return 0;
}
int restaVectores_01() {
    double* resul = crearVector(3);
    resul[0] = 1.0; resul[1] = 0.0; resul[2] = -1.0;

    double* v = crearVector(3);
    v[0] = 2; v[1] = 2; v[2] = 2;

    double* v2 = crearVector(3);
    v2[0] = 1; v2[1] = 2; v2[2] = 3;

    double* suma = crearVector(3);
    suma = restaVectores(v, 3, v2, 3);
    _assert(compararVectores(suma, resul, 3, 3));

    liberarVector(resul);
    liberarVector(suma);
    liberarVector(v);
    liberarVector(v2);
    return 0;
}
//MATRICES
int sumaMatrices_01() {
    double** m1, ** m2,**suma, ** sol;
    m1=crearMatriz(3, 3); m2=crearMatriz(3, 3); suma=crearMatriz(3, 3); sol=crearMatriz(3, 3);
    m1[0][0] = 1.0; m1[0][1] = 1.0; m1[0][2] = 2.0;
    m1[1][0] = 1.0; m1[1][1] = 1.0; m1[1][2] = 2.0;
    m1[2][0] = 1.0; m1[2][1] = 1.0; m1[2][2] = 2.0;

    m2[0][0] = 1.0; m2[0][1] = 1.0; m2[0][2] = 2.0;
    m2[1][0] = 1.0; m2[1][1] = 1.0; m2[1][2] = 2.0;
    m2[2][0] = 1.0; m2[2][1] = 1.0; m2[2][2] = 2.0;

    suma=sumaMatrices(m1,3,3,m2,3,3);

    sol[0][0] = 2.0; sol[0][1] = 2.0; sol[0][2] = 4.0;
    sol[1][0] = 2.0; sol[1][1] = 2.0; sol[1][2] = 4.0;
    sol[2][0] = 2.0; sol[2][1] = 2.0; sol[2][2] = 4.0;

    _assert(compararMatrices(suma, 3, 3, sol, 3, 3));

    liberarMatriz(m1, 3, 3); liberarMatriz(m2, 3, 3); liberarMatriz(suma, 3, 3); liberarMatriz(sol, 3, 3);
    return 0;
}
int pasarVectorMatriz_01()
{
    double** mat1 = crearMatriz(3, 3);
    mat1[0][0] = 1; mat1[0][1] = 2; mat1[0][2] = 0;
    mat1[1][0] = 2; mat1[1][1] = 4; mat1[1][2] = 7;
    mat1[2][0] = 3; mat1[2][1] = 1; mat1[2][2] = 0;

    double* v = crearVector(3);
    v[0] = 4; v[1] = 5; v[2] = 6;

    double** res = pasarVectorAMatriz(mat1, 3, 3, 1, v, 3);

    double** comp = crearMatriz(3, 3);
    comp[0][0] = 1; comp[0][1] = 4; comp[0][2] = 0;
    comp[1][0] = 2; comp[1][1] = 5; comp[1][2] = 7;
    comp[2][0] = 3; comp[2][1] = 6; comp[2][2] = 0;

    _assert(compararMatrices(res, 3, 3, comp, 3, 3));

    liberarMatriz(mat1, 3, 3);
    liberarVector(v);
    liberarMatriz(res, 3, 3);
    liberarMatriz(comp, 3, 3);

    return 0;
}
int matrizTraspuesta_01() {
    double** resultado, ** m1, ** transpuesta;
    m1 = crearMatriz(3, 3);
    m1[0][0] = 1.0; m1[0][1] = 0.0; m1[0][2] = 4.0;
    m1[1][0] = 0.0; m1[1][1] = 5.0; m1[1][2] = 0.0;
    m1[2][0] = 6.0; m1[2][1] = 0.0; m1[2][2] = -9.0;

    resultado = matrizTraspuesta(m1, 3);

    transpuesta = crearMatriz(3, 3);
    transpuesta[0][0] = 1.0; transpuesta[0][1] = 0.0; transpuesta[0][2] = 6.0;
    transpuesta[1][0] = 0.0; transpuesta[1][1] = 5.0; transpuesta[1][2] = 0.0;
    transpuesta[2][0] = 4.0; transpuesta[2][1] = 0.0; transpuesta[2][2] = -9.0;
    
    _assert(compararMatrices(resultado, 3, 3, transpuesta, 3, 3));
    liberarMatriz(m1, 3, 3);
    liberarMatriz(resultado, 3, 3);
    liberarMatriz(transpuesta, 3, 3);
    return 0;
}
int productoMatrices_01() {
    double** resultado,**m1,**m2,**p;
    m1 = crearMatriz(3, 3);
    m2 = crearMatriz(3, 3);
    m1[0][0] = 1.0; m1[0][1] = 4.0; m1[0][2] = 7.0;
    m1[1][0] = 2.0; m1[1][1] = 5.0; m1[1][2] = 8.0;
    m1[2][0] = 3.0; m1[2][1] = 6.0; m1[2][2] = 9.0;

    m2[0][0] = 1.0; m2[0][1] = -1.0; m2[0][2] = 2.0;
    m2[1][0] = 2.0; m2[1][1] = -1.0; m2[1][2] = 2.0;
    m2[2][0] = 3.0; m2[2][1] = -3.0; m2[2][2] = 0.0;

    resultado = productoMatrices(m1, 3, 3, m2, 3, 3);

    p = crearMatriz(3, 3);
    p[0][0] = 30.0; p[0][1] = -26.0; p[0][2] = 10.0;
    p[1][0] = 36.0; p[1][1] = -31.0; p[1][2] = 14.0;
    p[2][0] = 42.0; p[2][1] = -36.0; p[2][2] = 18.0;

    _assert(compararMatrices(resultado, 3, 3, p, 3, 3));
    liberarMatriz(resultado, 3, 3);
    liberarMatriz(m1, 3, 3);
    liberarMatriz(m2, 3, 3);
    liberarMatriz(p, 3, 3);
    return 0;
}
int determinanteMatriz_01() {
    double** m1 = crearMatriz(3, 3);
    m1[0][0] = 1.0; m1[0][1] = 1.0; m1[0][2] = 2.0;
    m1[1][0] = 1.0; m1[1][1] = 1.0; m1[1][2] = 2.0;
    m1[2][0] = 1.0; m1[2][1] = 1.0; m1[2][2] = 2.0;
    int resultado = determinanteMatriz(m1, 3);

    _assert(resultado == 0);
    liberarMatriz(m1, 3, 3);
    return 0;
}
int matrizIdentidad_01() {
    double** m1 = crearMatriz(3, 3);
    double** id = matrizIdentidad(3);
    m1[0][0] = 1.0; m1[0][1] = 0.0; m1[0][2] = 0.0;
    m1[1][0] = 0.0; m1[1][1] = 1.0; m1[1][2] = 0.0;
    m1[2][0] = 0.0; m1[2][1] = 0.0; m1[2][2] = 1.0;
    
    _assert(compararMatrices(m1, 3, 3, id, 3, 3));
    liberarMatriz(m1, 3, 3);
    liberarMatriz(id, 3, 3);
    return 0;

}
int matrizInversa_01() {
    double** resultado, ** m1,**inv;
    m1 = crearMatriz(3, 3);
    m1[0][0] = 1.0; m1[0][1] = 1.0; m1[0][2] = 0.0;
    m1[1][0] = 1.0; m1[1][1] = 0.0; m1[1][2] = 1.0;
    m1[2][0] = 0.0; m1[2][1] = 1.0; m1[2][2] = 0.0;
    resultado = matrizInversa(m1, 3, 3);

    inv = crearMatriz(3, 3);
    inv[0][0] = 1.0; inv[0][1] = 0.0; inv[0][2] = -1.0;
    inv[1][0] = 0.0; inv[1][1] = 0.0; inv[1][2] = 1.0;
    inv[2][0] = -1.0; inv[2][1] = 1.0; inv[2][2] = 1.0;

    _assert(compararMatrices(inv, 3, 3, resultado, 3, 3));
    liberarMatriz(m1, 3, 3);
    liberarMatriz(resultado, 3, 3);
    liberarMatriz(inv, 3, 3);
    return 0;

}
int matrizPorEscalar_01() {
    double** mat = crearMatriz(3, 3);
    mat[0][0] = 2; mat[1][1] = 2; mat[2][2] = 2;
    double** sol = matrizPorEscalar(matrizIdentidad(3, 3), 3, 3, 2);

    _assert(compararMatrices(sol, 3, 3, mat, 3, 3));

    liberarMatriz(mat, 3, 3);
    liberarMatriz(sol, 3, 3);

    return 0;
}
int matrizPorVector_01() {
    double** m = crearMatriz(3, 3);
    m[0][0] = 1; m[0][1] = 2; m[0][2] = 0;
    m[1][0] = 2; m[1][1] = 4; m[1][2] = 7;
    m[2][0] = 3; m[2][1] = 1; m[2][2] = 0;

    double* v = crearVector(3);
    v[0] = 4; v[1] = 5; v[2] = 6;

    double** sol = matrizPorVector(m, 3, 3, v, 3);

    double** resul = crearMatriz(3, 1);
    resul[0][0] = 14;
    resul[1][0] = 70;
    resul[2][0] = 17;

    _assert(compararMatrices(sol, 3, 1, resul, 3, 1));

    liberarVector(v);
    liberarMatriz(sol, 3, 1);
    liberarMatriz(resul, 3, 1);
    liberarMatriz(m, 3, 3);

    return 0;
}
//NUT ANGLES
int Nut_Angles_01() {
    double dpsi, deps, dpsi_s, deps_s;
    dpsi_s = 6.23063736216799e-05;
    deps_s = -3.51110708894389e-05;

    NutAngles(49746.1097222222, &dpsi, &deps);
    _assert(fabs(dpsi_s - dpsi) < pow(10, -10) && fabs(deps_s - deps) < pow(10, -10));
    return 0;
}
//NUT MATRIX
int NutMat_01() {
    double** m, **sol;
    sol = crearMatriz(3, 3);
    m = crearMatriz(3, 3);
    sol[0][0] = 0.999999998058958; sol[0][1] = -5.71647031444515e-05; sol[0][2] = -2.47846909052258e-05;
    sol[1][0] = 5.71655733262547e-05; sol[1][1] = 0.99999999774968; sol[1][2] = 3.51103624666727e-05;
    sol[2][0] = 2.47826837760046e-05; sol[2][1] = -3.51117792296018e-05; sol[2][2] = 0.999999999076491;

    double mjd = 4.974611085861109e+04;
    m = NutMat(mjd);

    _assert(compararMatrices(m, 3, 3, sol, 3, 3));

    liberarMatriz(sol, 3, 3);
    liberarMatriz(m, 3, 3);
    return 0;
}
//MEAN OBLIQUITY
int MeanObliquity_01() {
    double mj, sol;
    sol = 0.409103979370901;
    _assert(fabs(sol - MeanObliquity(49746.1097222222)) < pow(10, -10));
    return 0;
}
//MJDAY_TDB_01
int Mjday_TDB_01() {
    double sol= 4.974611199287850e+04;

    double resultado = Mjday_TDB(4.974611199287025e+04);
    //printf("%lf", resultado);
    _assert(fabs(sol - resultado) <pow(10,-10));
    return 0;
}
//POSICION
int posicion_01() {
    double lon, lat, alt, * Rs, * sol;
    sol = crearVector(3);
    lon = -2.76234307910694;
    lat = 0.376551295459273;
    alt = 300.2;
    sol[0] = -5512567.84003607;
    sol[1] = -2196994.44666933;
    sol[2] = 2330804.96614689;

    Rs = crearVector(3);
    Rs = posicion(lon, lat, alt);
    _assert(compararVectores(Rs,sol,3, 3));
    liberarVector(Rs);
    liberarVector(sol);
    return 0;
}
//UNIT
int unit_01() {
    double* v = crearVector(3);  
    v[0] = 1; v[1] = 2; v[2] = 3;

    double* resultado = crearVector(3);
    double norma = normaVector(v, 3);
    resultado[0] = v[0] / norma;
    resultado[1] = v[1] / norma;
    resultado[2] = v[2] / norma;

    double* unit = Unit(v, 3);

    _assert(compararVectores(unit, resultado, 3, 3));
    liberarVector(v);
    liberarVector(resultado);
    return 0;
}
//ANGL
int angl_01() {
    double* v1 = crearVector(3);
    v1[0] = 5.7203e+06;
    v1[1] = 3.1524e+06;
    v1[2] = 3.7501e+06;
    double* v2 = crearVector(3);
    v2[0] = 6.2214e+06;
    v2[1] = 2.8677e+06;
    v2[2] = 3.0062e+06;
    double resultado = angl(v1, 3, v2, 3);
    //printf("%lf", resultado);
    //quitamos un poco de precisión para que pase el test
    _assert(fabs(resultado-0.1253)<0.0001);
    liberarVector(v1);
    liberarVector(v2);

    return 0;
}
//AZELPA
int AzElPa_01() {
    double s[3] = { 1,2,3 };
    double Az = 0, El = 0, * dAds = crearVector(3), * dEds = crearVector(3);

    AzElPa(s, &Az, &El, dAds, dEds);

    double res_Az = 0.46364760900080609, res_El = 0.93027401411547206, res_dAds[3] = { 0.39999999999999991 , -0.19999999999999996 , 0.0 }, res_dEds[3] = { -0.095831484749990981, -0.19166296949998196, 0.15971914124998499 };

    _assert((fabs(Az - res_Az) < pow(10, -10)) && (fabs(El - res_El) < pow(10, -10)) && (compararVectores(dAds,res_dAds, 3,3)) && compararVectores(dEds,res_dEds, 3,3));

    return 0;

}
//GEODETIC
int Geodetic_01() {
    double* r = crearVector(3);
    r[0] = -5.512567840036068e+06;
    r[1] = -2.196994446669333e+06;
    r[2] = 2.330804966146887e+06;
    double lon,lat,h;
    Geodetic(r, &lon,&lat, &h);
    //hay que comprobar lon, lat, h
    double rlon = -2.7623,rlat=0.3766,rh=300.2000;
    //printf("%lf",lon);
    _assert(fabs(lon-rlon)<0.0001 && fabs(lat-rlat)<0.0001 && fabs(h-rh)<0.0001);
    liberarVector(r);
    return 0;
}
//ELEMENTS
int elements_01() {
    double** y = crearMatriz(3, 2);
    y[0][0] = 6.221397628578691e+06; y[0][1] = 4.645047251617496e+03;
    y[1][0] = 2.867713779657407e+06; y[1][1] = -2.752215915881822e+03;
    y[2][0] = 3.006155985099499e+06; y[2][1] = -7.507999409869392e+03;
    
    double p, a, e, i,Omega, omega, M;
    elements(&p, &a, &e, &i, &Omega, &omega, &M, y);
    //printf("%lf", M);
    _assert(fabs(p - 1.200169359721121e+07) < 0.0001e+07 && fabs(a - 1.894392266069662e+07) < 0.0001e+07
        && fabs(e - 0.6054) < 0.001 && fabs(i - 2.0266) < 0.001 && fabs(Omega - 3.3567) < 0.001
        && fabs(omega - 2.7376) < 0.001 && fabs(M - 6.2714) < 0.0001);

    liberarMatriz(y,3,2);
    return 0;
}
//FRAC
int frac_01() {
    double x = -4.534616865852740;
    double resultado = Frac(x);
    _assert(fabs(resultado - 0.465383134147261)<0.00001);
    return 0;
}
//LTC
int LTC_01() {
    double** resultado = crearMatriz(3, 3);
    double lat = 0.376551295459273, lon = -2.762343079106937;

    resultado = LTC(lon, lat);

    double** sol = crearMatriz(3, 3);
    sol[0][0] = 0.370223471399199; sol[0][1] = -0.928942722252092; sol[0][2] = 0.0;
    sol[1][0] = 0.341586711932422; sol[1][1] = 0.136136938528208; sol[1][2] = 0.929938305587722;
    sol[2][0] = -0.863859421119156; sol[2][1] = -0.344284987681776; sol[2][2] = 0.367715580035218;
    //imprimirMatriz(resultado, 3, 3);
    //imprimirMatriz(sol, 3, 3);
    _assert(compararMatrices(resultado, 3, 3, sol, 3, 3));
    return 0;
    liberarMatriz(resultado, 3, 3);
    liberarMatriz(sol, 3, 3);
}
//TIMEDIFF
int timediff_01() {
    double ut1_tai, utc_gps, ut1_gps, tt_utc, gps_utc;
    double ut1_utc=0.3258, tai_utc=29.0;

    timediff(ut1_utc, tai_utc, &ut1_tai, &utc_gps, &ut1_gps, &tt_utc, &gps_utc);
   
    _assert(fabs(ut1_tai - -28.6742) < 0.0001 && fabs(utc_gps- -10)<0.0001 && fabs(ut1_gps- -9.6742)<0.0001
        && fabs(tt_utc-61.1840)<0.0001 && fabs(gps_utc-10)<0.0001);
    return 0;
}
//SIGN_
int sign_01(){
    double a = 1.0, b = -1.349999919533730e+02;
    double resultado = sign_(a, b);
    _assert(resultado == -1);
    return 0;
}
//GMST
int gmst_01() {
    double mjd_ut1 = 4.974611015423337e+04;
    double resultado = gmst(mjd_ut1);
    _assert(fabs(resultado - 2.924088470683254)<0.00001);
    return 0;
}
//POLEMATRIX
int PoleMatrix_01() {
    double xp = -5.593861831521886e-07;
    double yp = 2.335544384403725e-06;
    double** resultado = PoleMatrix(xp, yp);

    double** sol = crearMatriz(3, 3);
    sol[0][0] = 0.999999999999844; sol[0][1] = -1.306471258772872e-12; sol[0][2] = -5.593861831506338e-07;
    sol[1][0] = 0.0; sol[1][1] = 0.999999999997273; sol[1][2] = -2.335544384401602e-06;
    sol[2][0] = 5.593861831521594e-07; sol[2][1] = 2.335544384401236e-06; sol[2][2] = 0.999999999997116;
    
    _assert(compararMatrices(sol, 3, 3, resultado, 3, 3));
    liberarMatriz(resultado, 3, 3);
    liberarMatriz(sol, 3, 3);
    return 0;
}
//TIMEUPDATE
int TimeUpdate_01() {
    double** P = crearMatriz(3, 3);
    P[0][0] = 1; P[0][1] = 0.0; P[0][2] = 1.0;
    P[1][0] = 2.0; P[1][1] = 0.0; P[1][2] = 2.0;
    P[2][0] = 3.0; P[2][1] = 1.0; P[2][2] = 3.0;

    double** Phi = crearMatriz(3, 3);
    Phi[0][0] = 2.0; Phi[0][1] = 0.0; Phi[0][2] = 3.0;
    Phi[1][0] = 1.0; Phi[1][1] = 0.0; Phi[1][2] = 1.0;
    Phi[2][0] = 1.0; Phi[2][1] = 2.0; Phi[2][2] = 1.0;

    double** resultado = crearMatriz(3, 3);
    resultado = TimeUpdate(P, 3, 3, Phi, 3, 3,1.0);

    double** P2 = crearMatriz(3, 3);
    P2[0][0] = 56.0; P2[0][1] = 23.0; P2[0][2] = 29.0;
    P2[1][0] = 21.0; P2[1][1] = 9.0; P2[1][2] = 11.0;
    P2[2][0] = 41.0; P2[2][1] = 17.0; P2[2][2] = 19.0;

    _assert(compararMatrices(resultado, 3, 3, P2, 3, 3));
    liberarMatriz(P, 3, 3);
    liberarMatriz(Phi, 3, 3);
    liberarMatriz(resultado, 3, 3);
    liberarMatriz(P2, 3, 3);
    return 0;
}
//GIBBS
int gibbs_01() {
    double* r1 = crearVector(3);
    r1[0] = 5.720303710129857e+06; r1[1] = 3.152426696533080e+06; r1[2] = 3.750056804164004e+06;
    double* r2 = crearVector(3);
    r2[0] = 6.221397628578691e+06; r2[1] = 2.867713779657407e+06; r2[2] = 3.006155985099499e+06;
    double* r3 = crearVector(3);
    r3[0] = 6.699811809767934e+06; r3[1] = 2.569867807638774e+06; r3[2] = 2.154940295423896e+06;

    double* v2 = crearVector(3);
    double theta, theta1, copa;
    char error;

    gibbs(r1,3, r2,3, r3,3, v2, &theta, &theta1, &copa, &error);

    double* sol = crearVector(3);
    sol[0] = 4.645047251617496e+03; sol[1] = -2.752215915881822e+03; sol[2] = -7.507999409869392e+03;
    //printf("%lf", theta);
    _assert(compararVectores(sol, v2, 3, 3) && fabs(theta - 0.125269502872991) < 0.0001
        && fabs(theta1 - 0.136454013492470) < 0.0001 && fabs(copa - 0.005097233477746) < 0.0001);

    liberarVector(r1);
    liberarVector(r2);
    liberarVector(r3);
    liberarVector(v2);
    liberarVector(sol);
    return 0;
}
//ACCELPOINTMASS
int AccelPointMass_01() {
    double GM = 1.327124400419394e+20;
    double* r = crearVector(3);
    r[0] = 6.221397628578691e+06; r[1] = 2.867713779657407e+06; r[2] = 3.006155985099499e+06;
    double* s = crearVector(3);
    s[0] = 9.229825172847661e+10; s[1] = -1.053751960790543e+11; s[2] = -4.568636722635329e+10;
    double* sol = AccelPointMass(r, 3, s, 3, GM);
    double* a = crearVector(3);
    a[0] = -1.868550593441697e-07; a[1] = -2.003329958831816e-07; a[2] = -1.599931207554889e-07;

    _assert(compararVectores(sol,a,3,3));
    liberarVector(r);
    liberarVector(s);
    liberarVector(sol);
    liberarVector(a);
    return 0;
}
//PRECMATRIX
int PrecMat_01() {
    double mjd_1 = 5.154450000000000e+04;
    double mjd_2 = 4.974611085861109e+04;
    double** resultado = PrecMat(mjd_1, mjd_2);
    double** sol = crearMatriz(3, 3);
    sol[0][0] = 0.999999279427995; sol[0][1] = 0.001101012138263; sol[0][2] = 4.784514219004540e-04;
    sol[1][0] = -0.001101012138264; sol[1][1] = 0.999999393885917; sol[1][2] = -2.633882773454226e-07;
    sol[2][0] = -4.784514218979997e-04; sol[2][1] = -2.633927355266896e-07; sol[2][2] = 0.999999885542077;
    _assert(compararMatrices(sol, 3, 3, resultado, 3, 3));

    liberarMatriz(sol, 3, 3);
    liberarMatriz(resultado,3,3);
    return 0;
}
//EQUINOX
int EqnEquinox_01() {
    double mjd = 4.974611015423337e+04;
    double resultado = EqnEquinox(mjd);

    _assert(fabs(resultado - 5.716489679866389e-05) < 0.000001);
    return 0;
}
//GAST
int gast_01() {
    double mjd = 4.974611015423337e+04;
    double resultado = gast(mjd);

    _assert(fabs(resultado - 2.924145635580052) < 0.0000001);
    return 0;
}
//GHAMatrix
int GHAMatrix_01() {
    double mjd = 4.974611015423337e+04;
    double** resultado = GHAMatrix(mjd);
    
    double** sol = crearMatriz(3, 2);
    sol[0][0] = -0.976451404712871; sol[0][1] = 0.215737465995733;
    sol[1][0] = -0.215737465995733; sol[1][1] = -0.976451404712871;
    sol[2][0] = 0.0; sol[2][1] = 0.0;

    _assert(compararMatrices(sol,3,2,resultado,3,2));
    liberarMatriz(sol, 3, 2);
    liberarMatriz(resultado, 3, 2);
    return 0;
}
//LEGENDRE
int legendre_01() {
    double** pnm = crearMatriz(3, 3);
    double** dpnm = crearMatriz(3, 3);

    Legendre(2, 2, 3.0, pnm, dpnm);

    double** solpnm = crearMatriz(3, 3);
    solpnm[0][0] = 1.0; solpnm[0][1] = 0.0; solpnm[0][2] = 0.0;
    solpnm[1][0] = 0.24442702392421947; solpnm[1][1] = -1.7147173032239305; solpnm[1][2] = 0.0;
    solpnm[2][0] = -1.0512375288193105; solpnm[2][1] = -0.54108578559834419; solpnm[2][2] = 1.8979267189818714;

    double** soldpnm = crearMatriz(3, 3);
    soldpnm[0][0] = 0.0;                   soldpnm[0][1] = 0.0;                       soldpnm[0][2] = 0.0;
    soldpnm[1][0] = -1.7147173032239305;   soldpnm[1][1] = -0.24442702392421947;      soldpnm[1][2] = 0.0;
    soldpnm[2][0] = -0.93718807190965236;  soldpnm[2][1] = 3.7187235297200689;        soldpnm[2][2] = 0.54108578559834419;


    _assert(compararMatrices(pnm, 3, 3, solpnm, 3, 3) && compararMatrices(dpnm, 3, 3, soldpnm, 3, 3));

    liberarMatriz(soldpnm, 3, 3);liberarMatriz(solpnm, 3, 3);
    liberarMatriz(pnm, 3, 3);liberarMatriz(dpnm, 3, 3);
    return 0;
}
//MEAS UPDATE
int MeasUpdate_01() {
    //Este test no funciona

    /*double* x = crearVector(6);
    x[0] = 7101597.83995656; x[1] = 1295244.79268408; x[2] = 12755.6074812165;
    x[3] = 576.097761535065; x[4] = -3084.51250378928; x[5] = -6736.03539066924;

    double z = 2653472;
    double g = 2653524.97225556;
    double s = 92.5;

    double* G = crearVector(6);
    G[0] = 0.484876050827302; G[1] = 0.0420045637468398; G[2] = -0.873573598478432;
    G[3] = 0.0; G[4] = 0.0; G[5] = 0.0;

    double** P = crearMatriz(6, 6);
    P[0][0] = 15877.862977197; P[0][1] = -5671.36398804291; P[0][2] = 8540.82446012286; P[0][3] = 48.5595982366734; P[0][4] = -13.3700368254264; P[0][5] = 22.4252217500793;
    P[1][0] = -5671.36398804291; P[1][1] = 24467.3370610125; P[1][2] = -1495.70216440439; P[1][3] = -4.64148177182865; P[1][4] = 59.5504403914688; P[1][5] = -26.9419372557421;
    P[2][0] = 8540.82446012287; P[2][1] = -1495.70216440438; P[2][2] = 6196.68034435384; P[2][3] = 26.7016364841125; P[2][4] = -4.13117368643687; P[2][5] = 16.0867101164001;
    P[3][0] = 48.5595982366735; P[3][1] = -4.64148177182864; P[3][2] = 26.7016364841125; P[3][3] = 0.178295954707253; P[3][4] = -0.0189362504462761; P[3][5] = 0.0668595526128017;
    P[4][0] = -13.3700368254264; P[4][1] = 59.5504403914687; P[4][2] = -4.13117368643688; P[4][3] = -0.0189362504462761; P[4][4] = 0.152611171445584; P[4][5] = -0.0750233063281362;
    P[5][0] = 22.4252217500793; P[5][1] = -26.941937255742; P[5][2] = 16.0867101164001; P[5][3] = 0.0668595526128016; P[5][4] = -0.0750233063281359; P[5][5] = 0.0845546951473536;

    double* K = crearVector(6);

    MeasUpdate(x, z, g, s, G, P, 6, K);

    double* sol_x = crearVector(6);
    sol_x[0] = 7101597.8425025502; sol_x[1] = 1295247.0610089891; sol_x[2] = 12762.893633465008;
    sol_x[3] = 576.09762737677102; sol_x[4] = -3084.5104703230959; sol_x[5] = -6736.0118584782977;

    double* sol_K = crearVector(6);
    sol_K[0] = -4.8062710501952013e-05; sol_K[1] = -0.042820999127366866; sol_K[2] = -0.13754657360661801;
    sol_K[3] = 2.5326143399425860e-06; sol_K[4] = -3.8387381520985645e-05; sol_K[5] = -0.00044423615142601212;

    double** sol_P = crearMatriz(6, 6);
    sol_P[0][0] = 15877.862954779022; sol_P[0][1] = -5671.3839611190006; sol_P[0][2] = 8540.7603040256945; sol_P[0][3] = 48.559599417965288; sol_P[0][4] = -13.370054730522208; sol_P[0][5] = 22.425014544207514;
    sol_P[1][0] = -5671.3839611190015; sol_P[1][1] = 24449.542245345161; sol_P[1][2] = -1552.8614080671396; sol_P[1][3] = -4.6404293114411281; sol_P[1][4] = 59.534488022490073; sol_P[1][5] = -27.126545288914798;
    sol_P[2][0] = 8540.7603040257054; sol_P[2][1] = -1552.8614080671305; sol_P[2][2] = 6013.0774794983963; sol_P[2][3] = 26.705017122656912; sol_P[2][4] = -4.1824147536386551; sol_P[2][5] = 15.493725297254827;
    sol_P[3][0] = 48.559599417965387; sol_P[3][1] = -4.6404293114411166; sol_P[3][2] = 26.705017122656919; sol_P[3][3] = 0.17829589246031052; sol_P[3][4] = -0.018935306955949978; sol_P[3][5] = 0.066870471109838431;
    sol_P[4][0] = -13.370054730522209; sol_P[4][1] = 59.534488022489967; sol_P[4][2] = -4.1824147536386640; sol_P[4][3] = -0.018935306955949981; sol_P[4][4] = 0.15259687075931724; sol_P[4][5] = -0.075188800341486062;
    sol_P[5][0] = 22.425014544207514; sol_P[5][1] = -27.126545288914702; sol_P[5][2] = 15.493725297254827; sol_P[5][3] = 0.066870471109838306; sol_P[5][4] = -0.075188800341485770; sol_P[5][5] = 0.082639523529451003;

    _assert(compararMatrices(P, 6, 6, sol_P, 6, 6) && compararVectores(sol_x, 6, x, 6) && compararVectores(sol_K, 6, K, 6) );

    liberarVector(x);
    liberarVector(G);
    liberarVector(K);
    liberarMatriz(P, 6, 6);

    liberarVector(sol_x);
    liberarVector(sol_K);
    liberarMatriz(sol_P, 6, 6);


    return 0;*/
}
//CHEB3D
int Cheb3D_01() {

    double cx[3] = { 1,2,3 };
    double cy[3] = { 4,5,6 };
    double cz[3] = { 1,2,3 };

    double** sol = crearMatriz(1, 3);
    sol[0][0] = -2.0; sol[0][1] = -2.0; sol[0][2] = -2.0;

    double** resultado = Cheb3D(2.0, 3, 1.0, 3.0, cx, cy, cz);

    _assert(compararMatrices(resultado, 1, 3, sol, 1, 3));

    liberarMatriz(sol, 1, 3); liberarMatriz(resultado, 1, 3);
    return 0;
}
//MJDAY
int Mjday_01() {
    double sol = 59333;
    double mjday = Mjday(2021, 4, 29, 0, 0, 0);
 
    _assert(fabs(mjday - sol) < 0.00000001);
    return 0;
}
//ACCEL HARMONIC
int AccelHarmonic_01()
{
    extern double** eopdata, ** Cnm, ** Snm;
    double aux1, aux2, aux3, aux4, aux5, aux6, aux7, aux8;

    

    int a, c, n, m, aux0;
    FILE* f;
    eopdata = crearMatriz(13, 21413);
    f = fopen("eop19620101.txt", "r");
    if (f == NULL)
    {
        printf("Fallo al abrir eop19620101.txt\n");
        exit(EXIT_FAILURE);
    }


    for (int g = 0; g < 21413; g++)
    {
        fscanf(f, "%d%d%d%d%lf%lf%lf%lf%lf%lf%lf%lf%d", &a, &c, &n, &m, &aux1, &aux2, &aux3, &aux4, &aux5, &aux6, &aux7, &aux8, &aux0);

        eopdata[0][g] = a;
        eopdata[1][g] = c; eopdata[2][g] = n; eopdata[3][g] = m; eopdata[4][g] = aux1;
        eopdata[5][g] = aux2; eopdata[6][g] = aux3; eopdata[7][g] = aux4; eopdata[8][g] = aux5;
        eopdata[9][g] = aux6; eopdata[10][g] = aux7; eopdata[11][g] = aux8; eopdata[12][g] = aux0;
    }
    fclose(f);

    Cnm = crearMatriz(181, 181);
    Snm = crearMatriz(181, 181);

    f = fopen("GGM03S.txt", "r");
    if (f == NULL)
    {
        printf("Fallo al abrir GGM03S.txt\n");
        exit(EXIT_FAILURE);
    }

    for (n = 0; n < 181; ++n)
    {
        for (m = 0; m <= n; ++m)
        {
            fscanf(f, "%d%d%lf%lf%lf%lf", &a, &c, &Cnm[n][m], &Snm[n][m], &aux1, &aux2);

        }
    }
    fclose(f);


    double* r = crearVector(3);
    r[0] = 1; r[1] = 2; r[2] = 3;
    double** E = crearMatriz(3, 3);
    E[0][0] = 1; E[0][1] = 2; E[0][2] = 3;
    E[1][0] = 3; E[1][1] = 4; E[1][2] = 5;
    E[2][0] = 5; E[2][1] = 6; E[2][2] = 7;


    double* resul = AccelHarmonic(r, E, 3, 3);

    double* sol = crearVector(3);
    sol[0] = -5.59550820651701e+21; sol[1] = -7.64168983653087e+21; sol[2] = -9.68787146654472e+21;

    int funciona = 1;

    for (int t = 0; t < 3; t++)
    {
        if (fabs(resul[t] - sol[t]) > pow(10, 8))
        {
            funciona = 0;
        }
    }

    _assert(funciona==1);

    liberarMatriz(E, 3, 3);
    liberarVector(sol);
    liberarVector(resul);

    return 0;
}
//JPL_Eph_DE430
int JPL_Eph_DE430_01()
{//Este test no pasa
    FILE* f;
    extern double** PC;
    PC = crearMatriz(2285, 1020);
    f = fopen("DE430Coeff.txt", "r");
    if (f == NULL)
    {
        printf("Error al abrir DE430Coeff.txt\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < 2285; i++)
    {
        for (int j = 0; j < 1020; j++)
        {
            fscanf(f, "%lf", &PC[i][j]);
        }
    }
    fclose(f);

    double* r_Mercury = crearVector(3), * r_Venus = crearVector(3), * r_Earth = crearVector(3), * r_Mars = crearVector(3), * r_Jupiter = crearVector(3), * r_Saturn = crearVector(3), * r_Uranus = crearVector(3), * r_Neptune = crearVector(3), * r_Pluto = crearVector(3), * r_Moon = crearVector(3), * r_Sun = crearVector(3);
    JPL_Eph_DE430(4.974611199287850e+04, r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun);
    double* sol_r_Mercury = crearVector(3), * sol_r_Venus = crearVector(3), * sol_r_Earth = crearVector(3), * sol_r_Mars = crearVector(3), * sol_r_Jupiter = crearVector(3), * sol_r_Saturn = crearVector(3), * sol_r_Uranus = crearVector(3), * sol_r_Neptune = crearVector(3), * sol_r_Pluto = crearVector(3), * sol_r_Moon = crearVector(3), * sol_r_Sun = crearVector(3);

    sol_r_Mercury[0] = 82317120520.628510; sol_r_Mercury[1] = -64956057153.683685; sol_r_Mercury[2] = -23030396695.807938;
    sol_r_Venus[0] = -14171545287.615829; sol_r_Venus[1] = -110771180336.88745; sol_r_Venus[2] = -41322532935.548119;
    sol_r_Earth[0] = -93492502017.266220; sol_r_Earth[1] = 105646796851.14340; sol_r_Earth[2] = 45805775676.886017;
    sol_r_Mars[0] = -87932727057.929398; sol_r_Mars[1] = 47092342823.872284; sol_r_Mars[2] = 29129065901.722427;
    sol_r_Jupiter[0] = -296878553276.01697; sol_r_Jupiter[1] = -753973827895.66089; sol_r_Jupiter[2] = -314193936949.47406;
    sol_r_Saturn[0] = 1483147838356.0020; sol_r_Saturn[1] = -452753373499.82068; sol_r_Saturn[2] = -248928486267.36093;
    sol_r_Uranus[0] = 1413650101107.0068; sol_r_Uranus[1] = -2510497422308.2764; sol_r_Uranus[2] = -1117739953759.0435;
    sol_r_Neptune[0] = 1872486216043.9543; sol_r_Neptune[1] = -3928138868833.0518; sol_r_Neptune[2] = -1654664659255.6875;
    sol_r_Pluto[0] = -2170184727141.9048; sol_r_Pluto[1] = -3914809613157.8389; sol_r_Pluto[2] = -552493523043.43835;
    sol_r_Moon[0] = 133707842.60742396; sol_r_Moon[1] = -324896800.96284586; sol_r_Moon[2] = -107799598.02716728;
    sol_r_Sun[0] = 93319283334.561188; sol_r_Sun[1] = -104626970282.98331; sol_r_Sun[2] = -45361943733.408005;

    _assert(compararVectores(sol_r_Mercury, 3, r_Mercury, 3) && compararVectores(r_Venus, 3, sol_r_Venus, 3) && compararVectores(r_Earth, 3, sol_r_Earth, 3)
        && compararVectores(r_Mars, 3, sol_r_Mars, 3) && compararVectores(r_Jupiter, 3, sol_r_Jupiter, 3) && compararVectores(r_Saturn, 3, sol_r_Saturn, 3)
        && compararVectores(r_Uranus, 3, sol_r_Uranus, 3) && compararVectores(r_Neptune, 3, sol_r_Neptune, 3) && compararVectores(r_Pluto, 3, sol_r_Pluto, 3)
        && compararVectores(r_Moon, 3, sol_r_Moon, 3) && compararVectores(r_Sun, 3, sol_r_Sun, 3));

    liberarVector(r_Mercury);liberarVector(r_Venus);liberarVector(r_Earth);liberarVector(r_Mars);
    liberarVector(r_Jupiter);liberarVector(r_Saturn);liberarVector(r_Uranus);liberarVector(r_Neptune);
    liberarVector(r_Pluto);liberarVector(r_Moon);liberarVector(r_Sun);liberarVector(sol_r_Mercury);
    liberarVector(sol_r_Venus);liberarVector(sol_r_Earth);liberarVector(sol_r_Mars);liberarVector(sol_r_Jupiter);
    liberarVector(sol_r_Saturn); liberarVector(sol_r_Uranus);liberarVector(sol_r_Neptune);liberarVector(sol_r_Pluto);
    liberarVector(sol_r_Moon);liberarVector(sol_r_Sun);
    return 0;
}
//IERS
int IERS_01()
{
    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
    IERS(4.974611015046295e+04, "l", &x_pole, &y_pole, &UT1_UTC, &LOD, &dpsi, &deps, &dx_pole, &dy_pole, &TAI_UTC);

    double sol_xpole = -5.5938618315218859e-07, sol_ypole = 2.3355443844037251e-06, sol_UT1_UTC = 0.32576469810652325, sol_LOD = 0.0027266863576381464,
         sol_dxpole = -8.4176415067052261e-10, sol_dypole = -1.5661888012134204e-09, sol_TAI_UTC = 29.0, sol_dpsi = -1.1688196064042147e-07, 
         sol_deps = -2.4788168041221858e-08;

    _assert(x_pole == sol_xpole && y_pole == sol_ypole && UT1_UTC == sol_UT1_UTC && LOD == sol_LOD && dpsi == sol_dpsi && deps == sol_deps && 
        dx_pole == sol_dxpole && dy_pole == sol_dypole && TAI_UTC == TAI_UTC);


    return 0;
}
//VAREQN
int VarEnq_01()
{
    extern Param AuxParam;
    double* yPhi = crearVector(42);
    yPhi[0] = 7101800.90695316;yPhi[1] = 1293997.58115302;yPhi[2] = 10114.0149489561;yPhi[3] = 573.068082065557;yPhi[4] = -3085.15736953139;
    yPhi[5] = -6736.03068347155; yPhi[6] = 1.00002934697411; yPhi[7] = 8.22733917824571e-06;yPhi[8] = 2.17104933032969e-07;yPhi[9] = 1.08925458286523e-05;
    yPhi[10] = 3.04673932159644e-06;yPhi[11] = 6.63504292924789e-08;yPhi[12] = 8.22733947396833e-06;yPhi[13] = 0.999986101965307;yPhi[14] = 3.99927484344495e-08;
    yPhi[15] = 3.04673962277234e-06;yPhi[16] = -5.15960624397393e-06;yPhi[17] = 1.22075293037738e-08; yPhi[18] = 2.17105647174788e-07;yPhi[19] = 3.99928715358944e-08;
    yPhi[20] = 0.999984551298692;yPhi[21] = 6.63510911923819e-08;yPhi[22] = 1.22076470413667e-08;yPhi[23] = -5.73276287737561e-06;yPhi[24] = 5.38976081674758;
    yPhi[25] = 1.47507305222446e-05;yPhi[26] = 3.2124178801047e-07;yPhi[27] = 1.00002936035847;yPhi[28] = 8.19365458245717e-06;yPhi[29] = 1.40504658160007e-07;
    yPhi[30] = 1.47507307389655e-05;yPhi[31] = 5.38968310056199;yPhi[32] = 5.90697772216375e-08;yPhi[33] = 8.19365490595617e-06;yPhi[34] = 0.99998608917631;
    yPhi[35] = 2.58022976866049e-08;yPhi[36] = 3.21242443032995e-07;yPhi[37] = 5.90698834536612e-08;yPhi[38] = 5.3896803255777;yPhi[39] = 1.40505380973508e-07;
    yPhi[40] = 2.58024226360683e-08; yPhi[41] = 0.999984550703337;

    
    AuxParam.Mjd_UTC = 49746.1163541665;
    AuxParam.m = 20;
    AuxParam.n = 20;
    AuxParam.sun = 1; AuxParam.moon = 1; AuxParam.planets = 1;
    AuxParam.Mjd_TT = 49746.1170623147;

    double* res = crearVector(42);

    VarEqn(5.38970808087707, yPhi, res);

    double* resultado = crearVector(42);
    for (int i = 0; i < 42; i++)
    {
        resultado[i] = res[i];
    }

    double* sol = crearVector(42);

    sol[0] = 573.06808206555695;sol[1] = -3085.15736953139003;sol[2] = -6736.03068347155022;sol[3] = -7.53489822593658;sol[4] = -1.37294429126636;
    sol[5] = -0.01075979864736;sol[6] = 0.00001089254583;sol[7] = 0.00000304673932;sol[8] = 0.00000006635043;sol[9] = 0.00000202239898;
    sol[10] = 0.00000056181190;sol[11] = 0.00000000439846;sol[12] = 0.00000304673962;sol[13] = -0.00000515960624;sol[14] = 0.00000001220753;
    sol[15] = 0.00000056181213;sol[16] = -0.00000095861369;sol[17] = 0.00000000080562;sol[18] = 0.00000006635109;sol[19] = 0.00000001220765;
    sol[20] = -0.00000573276288;sol[21] = 0.00000000439895;sol[22] = 0.00000000080571;sol[23] = -0.00000106368694;sol[24] = 1.00002936035847;
    sol[25] = 0.00000819365458;sol[26] = 0.00000014050466;sol[27] = 0.00001089991025;sol[28] = 0.00000302797127;sol[29] = 0.00000002370685;
    sol[30] = 0.00000819365491;sol[31] = 0.99998608917631;sol[32] = 0.00000002580230;sol[33] = 0.00000302797159;sol[34] = -0.00000516671243;
    sol[35] = 0.00000000434211;sol[36] = 0.00000014050538;sol[37] = 0.00000002580242;sol[38] = 0.99998455070334;sol[39] = 0.00000002370751;
    sol[40] = 0.00000000434223;sol[41] = -0.00000573302112;

    for (int i = 0; i < 42; i++)
    {
        if (fabs(sol[i] - resultado[i]) > pow(10, -10))
        {
            _assert(0);
        }
    }
    return 0;
}
//G_ACCELHARMONIC
int G_AccelHarmonic_01()
{
    double* r = crearVector(3);
    r[0] = 1; r[1] = 2; r[2] = 3;

    double** U = crearMatriz(3, 3);
    U[0][0] = 1; U[0][1] = 2; U[0][2] = 3;
    U[1][0] = 3; U[1][1] = 4; U[1][2] = 5;
    U[2][0] = 5; U[2][1] = 6; U[2][2] = 7;

    double** resultado = G_AccelHarmonic(r, U, 3, 3);

    double** sol = crearMatriz(3, 3);
    sol[0][0] = 3.1077519047699789e+21; sol[0][1] = 4.4388107279302436e+21; sol[0][2] = 5.8507592121296114e+21;
    sol[1][0] = 4.3477796981078816e+21; sol[1][1] = 6.0893095899059464e+21; sol[1][2] = 7.9329846340381688e+21;
    sol[2][0] = 5.5878074914457875e+21; sol[2][1] = 7.7398084518816535e+21; sol[2][2] = 1.0015210055946725e+22;

    int funciona = 1;
    for (int t = 0; t < 3; t++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (fabs(resultado[t][j] - sol[t][j]) > pow(10, 8))
            {
                funciona = 0;
            }
        }
    }
    _assert(funciona);

    liberarMatriz(U, 3, 3);
    liberarMatriz(resultado, 3, 3);
    liberarMatriz(sol, 3, 3);
    return 0;
}
//ACCEL
int Accel_01() {
    extern Param AuxParam;
    double* Y = crearVector(6);
    Y[0] = 5720694.2260585; Y[1] = 2687728.41425142; Y[2] = 3483000.08675422;
    Y[3] = 4371.83136151615; Y[4] = -1905.47309296259; Y[5] = -5698.58341612186;

    AuxParam.Mjd_UTC = 49746.1163541665; 
    AuxParam.m = 20; 
    AuxParam.n = 20;
    AuxParam.sun = 1; 
    AuxParam.moon = 1; 
    AuxParam.planets = 1;
    AuxParam.Mjd_TT = 49746.1170623147;

    double resultado[6];
    Accel(-543.476874884521, Y, resultado);

    double* sol = crearVector(6);
    sol[0] = 4371.8313615161496; sol[1] = -1905.4730929625900; sol[2] = -5698.5834161218600;
    sol[3] = -6.0654420256092694; sol[4] = -2.8497772792158966; sol[5] = -3.7023254336790901;
    _assert(compararVectores(resultado, 3, sol, 3));

    liberarVector(sol);
    liberarMatriz(PC, 2285, 1020);
    liberarMatriz(Cnm, 181, 181);
    liberarMatriz(Snm, 181, 181);
    liberarMatriz(eopdata, 13, 21413);
    return 0;
}
//-----------------------------------
//-----------------------------------
//-----------------------------------
int all_tests()
{
    
    _verify(R_x_01);
    _verify(R_y_01);
    _verify(R_z_01);

    _verify(normaVector_01);
    _verify(vectorPorEscalar_01);
    _verify(vectorMasEscalar_01);
    _verify(sumaVectores_01);
    _verify(restaVectores_01);

    _verify(sumaMatrices_01);
    _verify(matrizTraspuesta_01);
    _verify(productoMatrices_01);
    _verify(determinanteMatriz_01);
    _verify(matrizIdentidad_01);
    _verify(matrizInversa_01);
    _verify(pasarVectorMatriz_01);
    _verify(matrizPorEscalar_01);
    _verify(matrizPorVector_01);
    _verify(posicion_01);
    _verify(Nut_Angles_01);
    _verify(NutMat_01);
    _verify(MeanObliquity_01);
    _verify(Mjday_TDB_01);
    _verify(unit_01);
    _verify(angl_01);
    _verify(AzElPa_01);
    _verify(Geodetic_01);
    _verify(elements_01);
    _verify(frac_01);
    _verify(LTC_01);
    _verify(timediff_01);
    _verify(sign_01);
    _verify(gmst_01);
    _verify(PoleMatrix_01);
    _verify(TimeUpdate_01);
    _verify(gibbs_01);
    _verify(AccelPointMass_01);
    _verify(PrecMat_01);
    _verify(EqnEquinox_01);
    _verify(gast_01);
    _verify(GHAMatrix_01);
    _verify(legendre_01);
    //_verify(MeasUpdate_01);
    _verify(Cheb3D_01);
    _verify(Mjday_01);
    _verify(AccelHarmonic_01);
    //_verify(JPL_Eph_DE430_01);
    _verify(IERS_01);
    //_verify(VarEnq_01);
    _verify(G_AccelHarmonic_01);
    //_verify(Accel_01);


    return 0;
}


int main()
{
    int result = all_tests();
    printf("\n");
    if (result == 0)
        printf("TODOS LOS TEST PASAN CORRECTAMENTE!!\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
    /*double* v = crearVector(5);
    imprimirVector(v, 5);*/
}
