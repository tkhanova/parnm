#pragma once
#include "stdio.h"

#define FLOAT_TYPE       double
#define ILU_OK                 0
#define MATRIX_NOT_SQUARTE     1
#define NO_MATRIX_IN_FILE      2
#define MISMATCH_MATRIX_FORMAT 3
#define FIRST_VALUE_ZERO       4
#define CANT_OPEN_FILE         5
#define MATRIX_NOT_SYMMETRIC   6
#define EPSILON        0.0000001

//���� ������
#define UNDEFINE_TYPE   -1
#define SYMMETRIC        1
#define UPPER_TRIANGULAR 10 
#define LOWER_TRIANGULAR 11 

struct crsMatrix
{
  int N;  // ������ ������� (N x N)
  int NZ; // ���-�� ��������� ���������

  // ������ �������� (������ NZ)
  double* Value;
  // ������ ������� �������� (������ NZ)
  int* Col;
  // ������ �������� ����� (������ N + 1)
  int* RowIndex; 
  ~crsMatrix()
  {
	  if (Value!=NULL)
	  {
		  delete[] Value;
		  Value = NULL;
	  }
	  if (Col!=NULL)
	  {
		  delete[] Col;
		  Col = NULL;
	  }
	  if (RowIndex!=NULL)
	  {
		  delete[] RowIndex;
		  RowIndex = NULL;
	  }
  }
};
