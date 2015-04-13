#pragma once

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "type.h"


/**
 * API
 *   int ilu0(int n, double* a, int* col, int* row, 
 *      double* luval, int* uptr)
 *   ilu0 - ���������� �������
 * INPUT
 *   int    n   - ������ �������
 *   double * a - �� ������� ��������
 *   int  * col - ������� �������
 *   int  * row - ������� ������ �����
 * OUTPUT
 *   double * luval - �������� L � U ����������� ������
 *   int    * uptr  - ������� ������������ ���������
 *                    � ������� luval
 * RETURN
 *   ������������ ��� ������
 *   0    - ���������� ��������� �������
 *   -(n + 1) - ����� ������ ��� �� ��������� 0
 **/
int ilu0(int n, double * a, int * col, int * row, 
     double * luval, int * uptr);

