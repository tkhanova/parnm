#pragma once

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "type.h"


/**
 * API
 *   int ilu0(int n, double* a, int* col, int* row, 
 *      double* luval, int* uptr)
 *   ilu0 - разложение матрицы
 * INPUT
 *   int    n   - размер матрицы
 *   double * a - не нулевые элементы
 *   int  * col - индексы колонок
 *   int  * row - индексы начала строк
 * OUTPUT
 *   double * luval - значения L и U разложенных матриц
 *   int    * uptr  - индексы диагональных элементов
 *                    в массиве luval
 * RETURN
 *   возвращается код ошибки
 *   0    - разложение выполнено успешно
 *   -(n + 1) - номер строки где на диагонале 0
 **/
int ilu0(int n, double * a, int * col, int * row, 
     double * luval, int * uptr);

