#include <stdio.h>
#include "routines.h"
#include "readMTX.h"
#include "ilu0.h"
#include "timer.hpp"
#include "validation.h"
#include "conditionNumber.h"

int main(int argc, char **argv)
{
  // ������ ����������
  char *matrixName;
  ParseArgv(argc, argv, matrixName);

  // ���������� ����������� ����������
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

  // ������ ������� �� �����
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

  // ���� ������� ������������ � ������ ������ ������� �������������,
  // �� ��������� �� �� ������
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

  
  // �������� ��������� ������ ilu
  time->reset();
  time->start();
  InitializeMatrix(matA->N, matA->NZ, lu);
 
  memcpy(lu.Col, matA->Col, matA->NZ * sizeof(int));
  memcpy(lu.RowIndex, matA->RowIndex, (matA->N + 1) * sizeof(int));
  

  // �������� LU ����������
  error = ilu0(matA->N, matA->Value, matA->Col, matA->RowIndex,
  lu.Value, uptr);
  time->stop();

  printf("ILU factorization time: %f\n", time->getElapsed());

  if(error != ILU_OK)
  {
    printf("ilu factorization %d\n", error);
    return -1;
  }

  // �������� ������������ ILU

  // ���������� ������� �� L � �
  LUmatrixSeparation(lu, uptr, L, U);

  // ����������� �������
  ProductSparseMatrix(L, U, M);

#if DEBUG_LEVEL > 2
  // ���������� �������� ������������ ���������
  PrintMatrixInNormalView(matA);
  PrintMatrixInNormalView(&lu);
  PrintMatrixInNormalView(&L);
  PrintMatrixInNormalView(&U);
  PrintMatrixInNormalView(&M);
#endif

  // �������� ������������ ��������� ���������� �������
  if(!structValidation(*matA, M))
  {
    printf("invalid struct of matrix M \n");
    return -2;
  }

  // ���������� ������� ������� ��������
  diff = MatrixCompare(*matA, M);

  printf("distinction value of matrix %f \n", diff);

  if(matA->N < 4000)
  {
    //�������� ����� ��������������
    double condA, condMA;
    condA = getConditionNumber(*matA);
    printf("Condition number of matrix        A: %lf\n", condA);
    condMA = getConditionNumber(L, U, *matA);
    printf("Condition number of matrix M^-1 * A: %lf\n", condMA);
  }
  return 0;
}