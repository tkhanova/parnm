#include "validation.h"

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
void LUmatrixSeparation(crsMatrix ilu, int *uptr, crsMatrix &L, crsMatrix &U)
{
  int countL, countU;
  int i, j, s, f, k;
  double *val;
  int    *col;
  countU = 0;
  for(i = 0; i < ilu.N; i++)
  {
    countU += (ilu.RowIndex[i+1] - uptr[i]);
  }
  countL = ilu.NZ + ilu.N - countU;
  InitializeMatrix(ilu.N, countL, L);
  InitializeMatrix(ilu.N, countU, U);
  // fill L matrix
  k = 0;
  val = L.Value;
  col = L.Col;
  L.RowIndex[0] = k;
  for(i = 0; i < ilu.N; i++)
  {
    s = ilu.RowIndex[i];
    f = uptr[i];
    for(j = s; j < f; j++)
    {
      val[k] = ilu.Value[j];
      col[k] = ilu.Col[j];
      k++;
    }
    val[k] = 1.0;
    col[k] = i;
    k++;
    L.RowIndex[i + 1] = k;
  }
  // fill U matrix
  k = 0;
  val = U.Value;
  col = U.Col;
  U.RowIndex[0] = k;
  for(i = 0; i < ilu.N; i++)
  {
    s = uptr[i];
    f = ilu.RowIndex[i + 1];
    for(j = s; j < f; j++)
    {
      val[k] = ilu.Value[j];
      col[k] = ilu.Col[j];
      k++;
    }
    U.RowIndex[i + 1] = k;
  }
}

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
void ProductSparseMatrix(crsMatrix &A, crsMatrix &B, crsMatrix &C)
{
  int n = A.N;

  // �������� ��������� ��� ������ ������� MKL
  // ��������������� ������� A � B � �������
  int i, j;
  for (i = 0; i < A.NZ; i++)
    A.Col[i]++;
  for (i = 0; i < B.NZ; i++)
    B.Col[i]++;
  for (j = 0; j <= n; j++)
  {
    A.RowIndex[j]++;
    B.RowIndex[j]++;
  }

  // ������������ �������, ����������� C = op(A) * B
  char trans = 'N'; // ������� � ���, op(A) = A - �� ����� ��������������� A

// ������ ��������, �������� �� ��, ��� ����� ���������� ������
// request = 0: ������ ��� �������������� ������� �.�. �������� �������
// ���� �� �� �����, ������� ������ ���������� ��� �������� ����������,
// ����������:
// 1) �������� ������ ��� ������� �������� ����� ic: "���-�� �����+1" ���������;
// 2) ������� ������� � ���������� request = 1 - � ������� ic ����� �������� 
//                                                         ��������� �������
// 3) �������� ������ ��� �������� c � jc 
//    (���-�� ��������� = ic[���-�� �����]-1)
// 4) ������� ������� � ���������� request = 2
  int request;

// ��� ���� ������������� ������: ���� ����������� ���������, ����� �� 
// ������������� ������� A, B � C. � ��� ��������������, ��� ��� �������
// �����������, �������������, �������� ������� "No-No-Yes", �������
// ������������� ������ ��������, ����� ����� ����� �� 1 �� 7 ������������
  int sort = 8;

// ���������� ��������� ���������.
// ������������ ������ ���� request = 0
  int nzmax = -1;

// ��������� ����������
  int info;

// ������� ������ ��� ������� � ������� C
  C.RowIndex = new int[n + 1];
// ��������� ���������� ��������� ��������� � ������� C
  request = 1;
  C.Value = 0;
  C.Col = 0;
  mkl_dcsrmultcsr(&trans, &request, &sort, &n, &n, &n, A.Value, A.Col, 
                  A.RowIndex, B.Value, B.Col, B.RowIndex, C.Value, C.Col,
                  C.RowIndex, &nzmax, &info);
  int nzc = C.RowIndex[n] - 1;
  C.Value = new double[nzc];
  C.Col   = new int[nzc];
// ��������� C = A * B
  request = 2;
  mkl_dcsrmultcsr(&trans, &request, &sort, &n, &n, &n, A.Value, A.Col,
                  A.RowIndex, B.Value, B.Col, B.RowIndex, C.Value, C.Col,
                  C.RowIndex, &nzmax, &info);
  C.N = n;
  C.NZ = nzc;

  // �������� � ����������� ���� ������� A, B � �
  for (i = 0; i < A.NZ; i++)
    A.Col[i]--;
  for (i = 0; i < B.NZ; i++)
    B.Col[i]--;
  for (i = 0; i < C.NZ; i++)
    C.Col[i]--;
  for (j = 0; j <= n; j++)
  {
    A.RowIndex[j]--;
    B.RowIndex[j]--;
    C.RowIndex[j]--;
  }
}

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
bool structValidation   (crsMatrix &A,  crsMatrix &M)
{
  int i, j, fA, fM;
  i = 0;
  j = 0;
  int k;
  
  if(A.N != M.N)
  {
    return false;
  }

  if(M.NZ < A.NZ)
  {
    return false;
  }

  for(k = 0; k < A.N; k++)
  {
    i  = M.RowIndex[k];
    fM = M.RowIndex[k + 1];
    j  = A.RowIndex[k];
    fA = A.RowIndex[k + 1];
    while((i < fM) && (j<fA))
    {
      if(M.Col[i] != A.Col[j])
      {
        i++;
      } 
      else 
      {
        j++;
        i++;
      }
    }
    if((i == fM) && (j != fA))
    {
      return false;
    }
  }
  return true;
}

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
double MatrixCompare (crsMatrix &A,  crsMatrix &M)
{
  double diff = 0.0;
    int i, j, fA, fM;
  i = 0;
  j = 0;
  int k;
  
  if(A.N != M.N)
  {
    return -1.0;
  }

  if(M.NZ < A.NZ)
  {
    return -2.0;
  }

  for(k = 0; k < A.N; k++)
  {
    i  = M.RowIndex[k];
    fM = M.RowIndex[k + 1];
    j  = A.RowIndex[k];
    fA = A.RowIndex[k + 1];
    while((i < fM) && (j<fA))
    {
      if(M.Col[i] != A.Col[j])
      {
        i++;
      } 
      else 
      {
        if(diff < fabs(M.Value[i] - A.Value[j]))
        {
          diff =  fabs(M.Value[i] - A.Value[j]);
        }
        j++;
        i++;
      }
    }
  }
  return diff; 
}