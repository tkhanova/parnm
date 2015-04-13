#pragma once

#include "routines.h"
#include <string.h>
#include <mkl_lapack.h>
#include <mkl_spblas.h>

/**
 * API
 *   double * CRStoGeneral(crsMatrix A)
 *   ������� ������� �� ������� CRS � �������
 * INPUT
 *   crsMatrix A - ������� � ������� CRS
 * OUTPUT
 *   
 * RETURN
 *   ������� � ������� ����
 **/
double * CRStoGeneral(crsMatrix &A);

/**
 * API
 *   double getConditionNumber(crsMatrix A)
 *   ������ ����� ��������������
 * INPUT
 *   crsMatrix A - ������� � ������� CRS
 * OUTPUT
 *   
 * RETURN
 *   ����� �������������� 
 **/
double getConditionNumber(crsMatrix &A);

/**
 * API
 *   double* multIntMatMat(crsMatrix TrMat, char tr, double *Matrix)
 *   ��������� ������� �������� � ����������� �� �������
 * INPUT
 *   crsMatrix TrMat - ����������� �������, ������� ���������� 
 *                     �������� � �������� �� �������
 *   char tr         - ��� ����������� �������
 *   double *Matrix  - ���������� ������ �������
 * OUTPUT
 *   
 * RETURN
 *   ���������� ������� 
 **/
double getConditionNumber(crsMatrix &L, crsMatrix &U, crsMatrix &A);



