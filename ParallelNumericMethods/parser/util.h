#ifndef __UTIL_H__
#define __UTIL_H__

#include "type.h"
#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include <time.h>
#include <vector>


// ������� ���������� ������� � ������� CRS (3 �������, ���������� � ����)
// �������� ������ ��� ���� Value, Col � RowIndex
// ���������� ����� �������� mtx
void InitializeMatrix	(int N, int NZ, crsMatrix &mtx);

// ����������� ������ ��� ����� mtx
void FreeMatrix(crsMatrix &mtx);

// ������� ����� imtx � omtx, ������� ������ ��� ���� Value, Col � RowIndex
void CopyMatrix(const crsMatrix& imtx, crsMatrix &omtx);

#endif // __UTIL_H__