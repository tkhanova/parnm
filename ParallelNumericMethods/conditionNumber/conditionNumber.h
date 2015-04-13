#pragma once

#include "routines.h"
#include <string.h>
#include <mkl_lapack.h>
#include <mkl_spblas.h>

/**
 * API
 *   double * CRStoGeneral(crsMatrix A)
 *   перевод матрицы из формата CRS в плотный
 * INPUT
 *   crsMatrix A - матрица в формате CRS
 * OUTPUT
 *   
 * RETURN
 *   матрица в плотном виде
 **/
double * CRStoGeneral(crsMatrix A);

/**
 * API
 *   double getConditionNumber(crsMatrix A)
 *   оценка числа обусловлености
 * INPUT
 *   crsMatrix A - матрица в формате CRS
 * OUTPUT
 *   
 * RETURN
 *   число обусловлености 
 **/
double getConditionNumber(crsMatrix A);

/**
 * API
 *   double* multIntMatMat(crsMatrix TrMat, char tr, double *Matrix)
 *   умножение матрицы обратной к треугольной на плотную
 * INPUT
 *   crsMatrix TrMat - треугольная матрица, которую необходимо 
 *                     обратить и умножить на плотную
 *   char tr         - тип треугольной матрицы
 *   double *Matrix  - умножаемая справа матрица
 * OUTPUT
 *   
 * RETURN
 *   умноженная матрица 
 **/
double getConditionNumber(crsMatrix L, crsMatrix U, crsMatrix A);



