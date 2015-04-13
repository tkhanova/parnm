#ifndef __SPARSE_MATRIX_OPERATION_PARALLEL__
#define __SPARSE_MATRIX_OPERATION_PARALLEL__

#include "type.h"
#include "error.h"

/**
 * API
 *   int StructTranspose(int n, int* column, int* row, 
 *                       int* &tColumn, int* &tRow);
 *   ���������������� ��������� �������
 * INPUT
 *   int  n        - ������ �������
 *   int* column   - CRS �������� �������
 *   int* row  
 * OUTPUT
 *   int* &tColumn - CRS �������� �������
 *   int* &tRow  
 * RETURN
 * ������������ ��� ������
**/
int StructTranspose(int n, int* column, int* row, int* &tColumn, int* &tRow);

#endif