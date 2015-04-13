#ifndef __SPARSE_MATRIX_OPERATION_PARALLEL__
#define __SPARSE_MATRIX_OPERATION_PARALLEL__

#include "type.h"
#include "error.h"

/**
 * API
 *   int StructTranspose(int n, int* column, int* row, 
 *                       int* &tColumn, int* &tRow);
 *   Транспонирование структуры матрицы
 * INPUT
 *   int  n        - размер матрицы
 *   int* column   - CRS описание матрицы
 *   int* row  
 * OUTPUT
 *   int* &tColumn - CRS описание матрицы
 *   int* &tRow  
 * RETURN
 * возвращается код ошибки
**/
int StructTranspose(int n, int* column, int* row, int* &tColumn, int* &tRow);

#endif