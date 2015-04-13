#include <stdio.h>
#include "routines.h"
#include "readMTX.h"
#include "ilup.h"
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
  int countL;
  int countU;
  int n;
  int p;

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

  n = matA->N;
  uptr = new int [n];

  // создание структуру матриц ilu
  lu.Col      = NULL;
  lu.RowIndex = NULL;
  lu.Value    = NULL;

  for(p = 0; p < 4; p++)
  {
    printf("\n\n############### ILU(%d) ################\n", p);
    time->reset();
    time->start();
    symbolicILUp(p, matA->N, matA->Col, matA->RowIndex, 
                 lu.Col, lu.RowIndex, lu.Value,
                 uptr, countL, countU);
    lu.N = matA->N;
    lu.NZ = countL + countU;
    time->stop();
    if(error != ILU_OK)
    {
      printf("ilu symbolic factorization %d\n", error);
      return -1;
    }
    printf("Count elements: maxL NZ(%d) L NZ(%d) U NZ(%d) \n", (n * n - n)/2, countL + n, countU);
    printf("ILU symbolic factorization time: %f\n", time->getElapsed());

    time->reset();
    time->start();
    
    // неполное LU разложение
    error = numericalILUp(matA->N, matA->Value, matA->Col, matA->RowIndex,
      lu.Col, lu.RowIndex, uptr, lu.Value);
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
    // вычисление степени отличия значений
    diff = MatrixCompare(*matA, M);

    printf("distinction value of matrix %f \n", diff);
 
    if(matA->N < 4000)
    {
      double condA, condMA;
      condA = getConditionNumber(*matA);
      printf("Condition number of matrix        A: %lf\n", condA);
      condMA = getConditionNumber(L, U, *matA);
      printf("Condition number of matrix M^-1 * A: %lf\n", condMA);
    }
    #if DEBUG_LEVEL > 1
    // визуальная проверка корректности алгоритма
      PrintMatrixInNormalView(matA);
      PrintMatrixInNormalView(&lu);
    #if DEBUG_LEVEL > 2
      PrintMatrixInNormalView(&L);
      PrintMatrixInNormalView(&U);
    #endif
      PrintMatrixInNormalView(&M);
    #endif
  }

  return 0;
}