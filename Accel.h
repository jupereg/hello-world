#pragma once
#include "IERS.h"
#include "NutMat.h"
#include "matriz.h"
#include "vector.h"
#include "Mjday_TDB.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "JPL_Eph_DE430.h"
#include "iodkf.h"
#include "AccelPointMass.h"
#include "AccelHarmonic.h"
#include "globales.h"

void Accel(double x, double* Y, double* dY);