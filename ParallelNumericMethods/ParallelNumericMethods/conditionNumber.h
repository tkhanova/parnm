#pragma once

#include <string.h>
#include "crs_general.h"
#include <mkl_lapack.h>
#include <mkl_spblas.h>

/**
 * API
 *   double * CRStoGeneral(CRS A)
 *   ������� ������� �� ������� CRS � �������
 * INPUT
 *   CRS A - ������� � ������� CRS
 * OUTPUT
 *   
 * RETURN
 *   ������� � ������� ����
 **/
double * CRStoGeneral(CRS A);

/**
 * API
 *   double getConditionNumber(CRS A)
 *   ������ ����� ��������������
 * INPUT
 *   CRS A - ������� � ������� CRS
 * OUTPUT
 *   
 * RETURN
 *   ����� �������������� 
 **/
double getConditionNumber(CRS A);

/**
 * API
 *   double* multIntMatMat(CRS TrMat, char tr, double *Matrix)
 *   ��������� ������� �������� � ����������� �� �������
 * INPUT
 *   CRS TrMat - ����������� �������, ������� ���������� 
 *                     �������� � �������� �� �������
 *   char tr         - ��� ����������� �������
 *   double *Matrix  - ���������� ������ �������
 * OUTPUT
 *   
 * RETURN
 *   ���������� ������� 
 **/
double getConditionNumber(CRS L, CRS U, CRS A);



