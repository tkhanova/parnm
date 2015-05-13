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
 *   символьная фаза ILU(p)
 * INPUT
 *   int    n     - размер матрицы
 *   матрица A
 *   int  * col   - индексы колонок матрицы a
 *   int  * row   - индексы начала строк матрицы a
 *   портрет матрицы LU
 *   int  * &lucol - индексы колонок матрицы lu
 *   int  * &lurow - индексы начала строк матрицы lu
 *   double * &luval - значения матрицы
 *   int  * uptr  - индексы диагональных элементов
 *                  в массиве luval
 * OUTPUT
 *   double * luval - значения L и U разложенных матриц
 *   int &countL    - размер матрицы L
 *   int &countU    - размер матрицы U
 * RETURN
 *   возвращается код ошибки
 **/
int symbolicILUp(int p, int n, int * col, int * row, 
                 int * &lucol, int * &lurow, double * &luval,
                 int * uptr, int &countL, int &countU);


int symbolicILUpWithMultiplication(int p, crsMatrix& A, 
	crsMatrix* LU, 
	int * uptr);

/**
 * API
 *   int numericalILUp(int n, double * a, int * col, int * row, 
 *                     int * lucol, int * lurow, int * uptr,
 *                     double * luval);
 *   численная фаза ILU(p)
 * INPUT
 *   int    n     - размер матрицы
 *   double * a   - не нулевые элементы
 *   int  * col   - индексы колонок матрицы a
 *   int  * row   - индексы начала строк матрицы a
 *   int  * lucol - индексы колонок матрицы lu
 *   int  * lurow - индексы начала строк матрицы lu
 *   int  * uptr  - индексы диагональных элементов
 *                  в массиве luval
 * OUTPUT
 *   double * luval - значения L и U разложенных матриц
 * RETURN
 *   возвращается код ошибки
 *   0    - разложение выполнено успешно
 *   -(n + 1) - номер строки где на диагонале 0
 **/
int numericalILUp(int p, int n, double * a, int * col, int * row, 
                  int * lucol, int * lurow, int * uptr,
                  double * luval, int NZ);