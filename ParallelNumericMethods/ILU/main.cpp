#include <stdio.h>
#include "routines.h"
#include "readMTX.h"
#include "ilu0.h"
#include "timer.hpp"
#include "validation.h"
#include "conditionNumber.h"

int main(int argc, char **argv)
{
  // разбор параметров
  char *matrixName;
  ParseArgv(argc, argv, matrixName);

  // объявление необходимых переменных
  crsMatrix readA;
  crsMatrix *matA;
  crsMatrix lu;
  crsMatrix L;
  crsMatrix U;
  crsMatrix M;
  int *uptr;
  int typeOfMatrix;
  int error;
  double diff;
  Stopwatch *time = createStopwatch();

  // чтение матрицы из файла
  printf("read matrix (%s) \n", matrixName);
  time->start();
  error = ReadMatrixFromFile(matrixName, &(readA.N), &(readA.NZ),
  &(readA.Col), &(readA.RowIndex), &(readA.Value),
  &(typeOfMatrix));

  if(error != ILU_OK)
  {
  printf("error read matrix %d\n", error);
  return error;
  }

  // если матрица симметричная и задана только верхним треугольником,
  // то дополняем ее до полной
  if(typeOfMatrix == UPPER_TRIANGULAR)
  {
  matA = UpTriangleMatrixToFullSymmetricMatrix(&readA);
  } 
  else 
  {
  matA = &readA;
  }
  time->stop();
  printf("read matrix from file time: %f\n", time->getElapsed());

  uptr = new int [matA->N];

  
  // создание структуру матриц ilu
  time->reset();
  time->start();
  InitializeMatrix(matA->N, matA->NZ, lu);
 
  memcpy(lu.Col, matA->Col, matA->NZ * sizeof(int));
  memcpy(lu.RowIndex, matA->RowIndex, (matA->N + 1) * sizeof(int));
  

  // неполное LU разложение
  error = ilu0(matA->N, matA->Value, matA->Col, matA->RowIndex,
  lu.Value, uptr);
  time->stop();

  printf("ILU factorization time: %f\n", time->getElapsed());

  if(error != ILU_OK)
  {
    printf("ilu factorization %d\n", error);
    return -1;
  }

  // проверка корректности ILU

  // разделение матрицы на L и Г
  LUmatrixSeparation(lu, uptr, L, U);

  // перемножаем матрицы
  ProductSparseMatrix(L, U, M);

#if DEBUG_LEVEL > 2
  // визуальная проверка корректности алгоритма
  PrintMatrixInNormalView(matA);
  PrintMatrixInNormalView(&lu);
  PrintMatrixInNormalView(&L);
  PrintMatrixInNormalView(&U);
  PrintMatrixInNormalView(&M);
#endif

  // проверка корректности структуры полученной матрицы
  if(!structValidation(*matA, M))
  {
    printf("invalid struct of matrix M \n");
    return -2;
  }

  // вычисление степени отличия значений
  diff = MatrixCompare(*matA, M);

  printf("distinction value of matrix %f \n", diff);

  if(matA->N < 4000)
  {
    //проверка числа обусловлености
    double condA, condMA;
    condA = getConditionNumber(*matA);
    printf("Condition number of matrix        A: %lf\n", condA);
    condMA = getConditionNumber(L, U, *matA);
    printf("Condition number of matrix M^-1 * A: %lf\n", condMA);
  }
  return 0;
}