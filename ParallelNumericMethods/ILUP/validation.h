#pragma once
#include "util.h"
#include "mkl_types.h"
#include "mkl_spblas.h"

/**
 * �������� ������������ ilu ����������
 */

/**
 * API
 *  void LUmatrixSeparation (crsMatrix ilu, int *uptr, 
 *                           crsMatrix &L, crsMatrix &U);
 *  ���������� ������� �� L � U �������
 * INPUT
 *   crsMatrix ilu  - ������� ilu � ����� ���������
 *   int    * uptr  - ������� ������������ ���������
 *                    � ������� �������� ilu
 * OUTPUT
 *   crsMatrix &L   - ���������� ������� L
 *   crsMatrix &U   - ���������� ������� U
 * RETURN
 */
void LUmatrixSeparation(const crsMatrix &ilu, int *uptr, crsMatrix &L, crsMatrix &U);
/**
 * API
 *  void ProductSparseMatrix (crsMatrix &A,  crsMatrix &B, 
 *                            crsMatrix &C);
 *  ��������� ����������� ������ C = A * B
 * INPUT
 *   crsMatrix &A   - ������� A
 *   crsMatrix &B   - ������� B
 * OUTPUT
 *   crsMatrix &C   - C = A * B
 * RETURN
 */
void ProductSparseMatrix (crsMatrix &A,  crsMatrix &B, 
                          crsMatrix &C);
/**
 * API
 *  bool structValidation (crsMatrix &A,  crsMatrix &M); 
 *  �������� ������������ ��������� �������������������
 * INPUT
 *   crsMatrix &A   - �������� �������
 *   crsMatrix &M   - �������������������
 * OUTPUT
 *   
 * RETURN
 *  ��������� �� ���������
 */
bool structValidation (crsMatrix &A,  crsMatrix &M);
/**
 * API
 *  double MatrixCompare (crsMatrix &A,  crsMatrix &M);
 *  ������� ������� ������� ������
 * INPUT
 *   crsMatrix &A   - �������� �������
 *   crsMatrix &M   - �������������������
 * OUTPUT
 *   
 * RETURN
 *  ������� �������
 */
double MatrixCompare (crsMatrix &A,  crsMatrix &M);
