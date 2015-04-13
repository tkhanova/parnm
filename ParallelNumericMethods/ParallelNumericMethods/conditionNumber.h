#pragma once

#include <string.h>
#include "crs_general.h"
#include <mkl_lapack.h>
#include <mkl_spblas.h>

/**
 * API
 *   double * CRStoGeneral(CRS A)
 *   перевод матрицы из формата CRS в плотный
 * INPUT
 *   CRS A - матрица в формате CRS
 * OUTPUT
 *   
 * RETURN
 *   матрица в плотном виде
 **/
double * CRStoGeneral(CRS A);

/**
 * API
 *   double getConditionNumber(CRS A)
 *   оценка числа обусловлености
 * INPUT
 *   CRS A - матрица в формате CRS
 * OUTPUT
 *   
 * RETURN
 *   число обусловлености 
 **/
double getConditionNumber(CRS A);

/**
 * API
 *   double* multIntMatMat(CRS TrMat, char tr, double *Matrix)
 *   умножение матрицы обратной к треугольной на плотную
 * INPUT
 *   CRS TrMat - треугольная матрица, которую необходимо 
 *                     обратить и умножить на плотную
 *   char tr         - тип треугольной матрицы
 *   double *Matrix  - умножаемая справа матрица
 * OUTPUT
 *   
 * RETURN
 *   умноженная матрица 
 **/
double getConditionNumber(CRS L, CRS U, CRS A);



