#pragma once
#include "util.h"
#include "mkl_types.h"
#include "mkl_spblas.h"

/**
 * Проверка корректности ilu разложения
 */

/**
 * API
 *  void LUmatrixSeparation (crsMatrix ilu, int *uptr, 
 *                           crsMatrix &L, crsMatrix &U);
 *  разделение матрицы на L и U матрицы
 * INPUT
 *   crsMatrix ilu  - матрицы ilu в одной структуре
 *   int    * uptr  - индексы диагональных элементов
 *                    в массиве значений ilu
 * OUTPUT
 *   crsMatrix &L   - отделенная матрица L
 *   crsMatrix &U   - отделенная матрица U
 * RETURN
 */
void LUmatrixSeparation (crsMatrix ilu, int *uptr, 
                         crsMatrix &L, crsMatrix &U);
/**
 * API
 *  void ProductSparseMatrix (crsMatrix &A,  crsMatrix &B, 
 *                            crsMatrix &C);
 *  умножение разреженных матриц C = A * B
 * INPUT
 *   crsMatrix &A   - матрица A
 *   crsMatrix &B   - матрица B
 * OUTPUT
 *   crsMatrix &C   - C = A * B
 * RETURN
 */
void ProductSparseMatrix (crsMatrix &A,  crsMatrix &B, 
                          crsMatrix &C);
/**
 * API
 *  bool structValidation (crsMatrix &A,  crsMatrix &M); 
 *  проверка корректности структуры предобуславлевотеля
 * INPUT
 *   crsMatrix &A   - исходная матрица
 *   crsMatrix &M   - предобуславлевотель
 * OUTPUT
 *   
 * RETURN
 *  правильна ли структура
 */
bool structValidation (crsMatrix &A,  crsMatrix &M);
/**
 * API
 *  double MatrixCompare (crsMatrix &A,  crsMatrix &M);
 *  подсчет степени отличия матриц
 * INPUT
 *   crsMatrix &A   - исходная матрица
 *   crsMatrix &M   - предобуславлевотель
 * OUTPUT
 *   
 * RETURN
 *  степень отличия
 */
double MatrixCompare (crsMatrix &A,  crsMatrix &M);
