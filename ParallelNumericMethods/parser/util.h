#ifndef __UTIL_H__
#define __UTIL_H__

#include "type.h"
#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include <time.h>
#include <vector>


// Создает квадратную матрицу в формате CRS (3 массива, индексация с нуля)
// Выделяет память под поля Value, Col и RowIndex
// Возвращает через параметр mtx
void InitializeMatrix	(int N, int NZ, crsMatrix &mtx);

// Освобождает память для полей mtx
void FreeMatrix(crsMatrix &mtx);

// Создает копию imtx в omtx, выделяя память под поля Value, Col и RowIndex
void CopyMatrix(crsMatrix imtx, crsMatrix &omtx);

#endif // __UTIL_H__