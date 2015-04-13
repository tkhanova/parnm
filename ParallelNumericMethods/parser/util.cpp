#include "util.h"

// Создает квадратную матрицу в формате CRS (3 массива, индексация с нуля)
// Выделяет память под поля Value, Col и RowIndex
// Возвращает через параметр mtx
void InitializeMatrix(int N, int NZ, crsMatrix &mtx)
{
  mtx.N = N;
  mtx.NZ = NZ;
  mtx.Value  = new double[NZ];
  mtx.Col    = new int[NZ];
  mtx.RowIndex = new int[N + 1];
}

// Освобождает память, выделенную под поля mtx
void FreeMatrix(crsMatrix &mtx)
{
  delete [] mtx.Value;
  delete [] mtx.Col;
  delete [] mtx.RowIndex;
  mtx.Value = NULL;
  mtx.Col = NULL;
  mtx.RowIndex = NULL;
}

// Создает копию imtx в omtx, выделяя память под поля Value, Col и RowIndex
void CopyMatrix(crsMatrix imtx, crsMatrix &omtx)
{
  // Инициализация результирующей матрицы
  int N  = imtx.N;
  int NZ = imtx.NZ;
  InitializeMatrix(N, NZ, omtx);
  // Копирование
  memcpy(omtx.Value   , imtx.Value   , NZ * sizeof(double));
  memcpy(omtx.Col   , imtx.Col   , NZ * sizeof(int));
  memcpy(omtx.RowIndex, imtx.RowIndex, (N + 1) * sizeof(int));
}