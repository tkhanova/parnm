#pragma once

#include <stdio.h>
#include <string.h>
#include <queue>
#include <math.h>
#include "type.h"
#include "sparseMatrixOperation.h"
/**
 * API
 *   int symbolicILUp(int p, int n, int * col, int * row, 
 *                    int * lucol, int * lurow, double * &luval,
 *                    int * uptr, int &countL, int &countU);
 *   ���������� ���� ILU(p)
 * INPUT
 *   int    n     - ������ �������
 *   ������� A
 *   int  * col   - ������� ������� ������� a
 *   int  * row   - ������� ������ ����� ������� a
 *   ������� ������� LU
 *   int  * &lucol - ������� ������� ������� lu
 *   int  * &lurow - ������� ������ ����� ������� lu
 *   double * &luval - �������� �������
 *   int  * uptr  - ������� ������������ ���������
 *                  � ������� luval
 * OUTPUT
 *   double * luval - �������� L � U ����������� ������
 *   int &countL    - ������ ������� L
 *   int &countU    - ������ ������� U
 * RETURN
 *   ������������ ��� ������
 **/
int symbolicILUp(int p, int n, int * col, int * row, 
                 int * &lucol, int * &lurow, double * &luval,
                 int * uptr, int &countL, int &countU);
/**
 * API
 *   int numericalILUp(int n, double * a, int * col, int * row, 
 *                     int * lucol, int * lurow, int * uptr,
 *                     double * luval);
 *   ��������� ���� ILU(p)
 * INPUT
 *   int    n     - ������ �������
 *   double * a   - �� ������� ��������
 *   int  * col   - ������� ������� ������� a
 *   int  * row   - ������� ������ ����� ������� a
 *   int  * lucol - ������� ������� ������� lu
 *   int  * lurow - ������� ������ ����� ������� lu
 *   int  * uptr  - ������� ������������ ���������
 *                  � ������� luval
 * OUTPUT
 *   double * luval - �������� L � U ����������� ������
 * RETURN
 *   ������������ ��� ������
 *   0    - ���������� ��������� �������
 *   -(n + 1) - ����� ������ ��� �� ��������� 0
 **/
int numericalILUp(int n, double * a, int * col, int * row, 
                  int * lucol, int * lurow, int * uptr,
                  double * luval);