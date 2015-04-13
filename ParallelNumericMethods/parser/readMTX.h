#pragma once
#include <stdio.h>
#include "routines.h"

/**
* API
*   int ReadMTXFile(FILE* matrixFile, int* n, int* countNotZero, int** columns, 
*          int** rows, FLOAT_TYPE** values, int** countValueInRow,
*          int* isOrdered, int* isStoreOnlyUpperTriangle, 
*          int* isStoreOnlyLowTriangle)
*   ������ ������� �� mtx ������� ����� � ������������ ������ �������� ������
* INPUT
*   FILE* matrixFile   - ���� �������
* OUTPUT
*   int  n         - ������ �������
*   int* countNotZero  - ���-�� ��������� ���������
*   int** columns    - ������� � ������������ �������
*   int** rows
*   FLOAT_TYPE** values
*   int** countValueInRow - ���������� �������� � ������
*   int* isOrdered    - ����������� �� �������
*   int* isStoreOnlyUpperTriangle - �������� �� ������ ������� �����������
*   int* isStoreOnlyLowTriangle   - �������� �� ������ ������ �����������
* RETURN
*   ������������ ��� ������
**/
int ReadMTXFile(FILE* matrixFile, int* n, int* countNotZero, int** columns, 
        int** rows, FLOAT_TYPE** values, int** countValueInRow,
        int* isOrdered, int* isStoreOnlyUpperTriangle, 
        int* isStoreOnlyLowTriangle);

/**
* API
*   CheckSymmetric(int countNotZero, int* columns, int* rows, 
*          FLOAT_TYPE* values)
*   �������� ������� �� ��������������
* INPUT
*   int countNotZero   - ����� ��������� ��������� � ������� 
* OUTPUT
*   int* columns     - ������� � ������������ �������
*   int* rows
*   FLOAT_TYPE* values
* RETURN
*   ������������ ����������� ������� ��� ���
**/
int CheckSymmetric(int countNotZero, int* columns, int* rows, 
           FLOAT_TYPE* values);
/**
* API
*   SortMatrix(int n, int* column, int* row, FLOAT_TYPE* val);
*   �������������� �������
* INPUT
*   int  n         - ������ �������
*   int* column      - CRS �������� �������
*   int* row  
*   FLOAT_TYPE* val
* OUTPUT
*   int* column      - CRS �������� �������
*   int* row  
*   FLOAT_TYPE* val
* RETURN
* ������������ ����� ������
**/
void SortMatrix(int n, int* column, int* row, FLOAT_TYPE* val);

/**
* API
*   ReadMatrixFromFile(char* matrixName, int* n, int** column, 
*            int** row, FLOAT_TYPE** val)
*   ������ ������� �� mtx ������� ����� � crs ������ �������� ������
* INPUT
*   char* matrixName
* OUTPUT
*   int  n         - ������ �������
*   int** column     - CRS �������� �������
*   int** row
*   FLOAT_TYPE** val
* RETURN
*   ������������ ��� ������
**/
int ReadMatrixFromFile(char* matrixName, int* n, int* nz, int** column, 
             int** row, FLOAT_TYPE** val, int *typeOfMatrix);

