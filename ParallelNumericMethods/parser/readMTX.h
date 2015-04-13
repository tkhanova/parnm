#pragma once
#include <stdio.h>
#include "routines.h"

/**
* API
*   int ReadMTXFile(FILE* matrixFile, int* n, int* countNotZero, int** columns, 
*          int** rows, FLOAT_TYPE** values, int** countValueInRow,
*          int* isOrdered, int* isStoreOnlyUpperTriangle, 
*          int* isStoreOnlyLowTriangle)
*   чтение матрицы из mtx формата файла в координатный формат хранения матриц
* INPUT
*   FILE* matrixFile   - файл матрицы
* OUTPUT
*   int  n         - размер матрицы
*   int* countNotZero  - кол-во ненулевых элементов
*   int** columns    - матрица в координатном формате
*   int** rows
*   FLOAT_TYPE** values
*   int** countValueInRow - количество значений в строке
*   int* isOrdered    - упорядочена ли матрица
*   int* isStoreOnlyUpperTriangle - хранится ли только верхний треугольник
*   int* isStoreOnlyLowTriangle   - хранится ли только нижний треугольник
* RETURN
*   возвращается код ошибки
**/
int ReadMTXFile(FILE* matrixFile, int* n, int* countNotZero, int** columns, 
        int** rows, FLOAT_TYPE** values, int** countValueInRow,
        int* isOrdered, int* isStoreOnlyUpperTriangle, 
        int* isStoreOnlyLowTriangle);

/**
* API
*   CheckSymmetric(int countNotZero, int* columns, int* rows, 
*          FLOAT_TYPE* values)
*   проверка матрицы на симметричность
* INPUT
*   int countNotZero   - число ненулевых элементов в матрице 
* OUTPUT
*   int* columns     - матрица в координатном формате
*   int* rows
*   FLOAT_TYPE* values
* RETURN
*   возвращается симметрична матрица или нет
**/
int CheckSymmetric(int countNotZero, int* columns, int* rows, 
           FLOAT_TYPE* values);
/**
* API
*   SortMatrix(int n, int* column, int* row, FLOAT_TYPE* val);
*   Упорядочивание матрицы
* INPUT
*   int  n         - размер матрицы
*   int* column      - CRS описание матрицы
*   int* row  
*   FLOAT_TYPE* val
* OUTPUT
*   int* column      - CRS описание матрицы
*   int* row  
*   FLOAT_TYPE* val
* RETURN
* возвращается время работы
**/
void SortMatrix(int n, int* column, int* row, FLOAT_TYPE* val);

/**
* API
*   ReadMatrixFromFile(char* matrixName, int* n, int** column, 
*            int** row, FLOAT_TYPE** val)
*   чтение матрицы из mtx формата файла в crs формат хранения матриц
* INPUT
*   char* matrixName
* OUTPUT
*   int  n         - размер матрицы
*   int** column     - CRS описание матрицы
*   int** row
*   FLOAT_TYPE** val
* RETURN
*   возвращается код ошибки
**/
int ReadMatrixFromFile(char* matrixName, int* n, int* nz, int** column, 
             int** row, FLOAT_TYPE** val, int *typeOfMatrix);

