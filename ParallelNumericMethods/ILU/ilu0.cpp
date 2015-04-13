#include "ilu0.h"

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
     double * luval, int * uptr)
{
  int j1, j2;     // граница текущей строки
  int jrow;       // номер текущего столбца
  int k, j, jj;   // счетчики циклов
  int *iw = NULL;
  int jw;
  double t1;

  iw = new int[n];
  memcpy(luval, a, row[n] * sizeof(double));

  memset(iw, 0, n * sizeof(int));
  
  for(k = 0; k < n; k++)
  {
  j1 = row[k];
  j2 = row[k + 1];
  for(j = j1; j < j2; j++)
  {
    iw[col[j]] = j;
  }
  for(j = j1; (j < j2) && (col[j] < k); j++)
  {
    jrow = col[j];
    t1 = luval[j] / luval[uptr[jrow]];
    luval[j] = t1;
    for(jj = uptr[jrow]+1; jj < row[jrow + 1]; jj++)
    {
    jw = iw[col[jj]];
    if(jw != 0)
    {
      luval[jw] = luval[jw] - t1 * luval[jj];
    }
    }
  }
  jrow = col[j];
  uptr[k] = j;
  if((jrow != k) || (fabs(luval[j]) < EPSILON))
  {
    break;
  }
  for(j = j1; j < j2; j++)
  {
    iw[col[j]] = 0;
  }
  }

  delete [] iw;
  if(k < n)
  return -(k+1);
  return 0;
}