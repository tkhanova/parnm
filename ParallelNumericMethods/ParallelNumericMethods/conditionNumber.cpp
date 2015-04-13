#include "conditionNumber.h"
#include "crs_general.h"
#include <math.h>

/**
 * API
 *   double MatrixNormO(double* Matrix, int size)
 *   ������ ����� �������
 * INPUT
 *   double* Matrix - ������� � ������� ����
 *   int     size   - ������ �������
 * OUTPUT
 *   
 * RETURN
 *   ������ ����� �������
 **/
double MatrixNormO(double* Matrix, int size)
{
  double norm = 0.0;
  double sum;
  for (int j = 0; j < size; j++)
  {
    sum = 0.0;
    for (int i = 0; i < size; i++) sum += fabs(Matrix[i*size + j]);
    if (sum > norm) norm = sum;
  }
  return norm;
}

/**
 * API
 *   double MatrixNormI(double* Matrix, int size)
 *   ����������� ����� �������
 * INPUT
 *   double* Matrix - ������� � ������� ����
 *   int     size   - ������ �������
 * OUTPUT
 *   
 * RETURN
 *   ����������� ����� �������
 **/
double MatrixNormI(double* Matrix, int size)
{
  double norm = 0.0;
  double sum;
  for (int i = 0; i < size; i++)
  {
    sum = 0.0;
    for (int j = 0; j < size; j++) sum += fabs(Matrix[i*size + j]);
    if (sum > norm) norm = sum;
  }
  return norm;
}

/**
 * API
 *   double * CRStoGeneral(CRS A)
 *   ������� ������� �� ������� CRS � �������
 * INPUT
 *   CRS A - ������� � ������� CRS
 * OUTPUT
 *   
 * RETURN
 *   ������� � ������� ����
 **/
double * CRStoGeneral(CRS& A)
{
  int i, j, s, f;
  double * mat;
  mat = new double[A.N * A.N];
  memset(mat, 0, A.N * A.N * sizeof(double));
  for(i = 0; i < A.N; i++)
  {
    s = A.row_ptr[i];
    f = A.row_ptr[i + 1];
    for(j = s; j < f; j++)
    {
      mat[i * A.N + A.col_ind[j]] = A.val[j];
    }
  }
  return mat;
}

/**
 * API
 *   double getConditionNumber(CRS A)
 *   ������ ����� ��������������
 * INPUT
 *   CRS A - ������� � ������� CRS
 * OUTPUT
 *   
 * RETURN
 *   ����� �������������� 
 **/
double getConditionNumber(CRS A)
{
  // ��������������� ����������
  double *Matrix;
  int* Size = &(A.N);
  int* ipiv = new int [A.N];
  double* ANorm; 
  double* work = new double [A.N*4];
  int* iwork = new int [A.N];
  int lda = A.N;
  int* Lda = &lda;
  int info = 0;
  int *Info = &info;
  double rcond;
  double* Rcond = &rcond;
  // ����� � ����� ���������������
  double norm, cond = -1.0;

  //�������� ������� � �������� ����
  Matrix = CRStoGeneral(A);

  //��������� �����
  norm = MatrixNormI(Matrix, A.N);

  //�LU ����������
  dgetrf(Size,Size,Matrix,Lda,ipiv,Info);
  
  //����� ����� ������ �����
  ANorm = &norm;
  dgecon("O",Size,Matrix,Lda,ANorm,Rcond,work,iwork,Info);
  cond = 1.0/(rcond);

  delete [] Matrix;
  delete [] ipiv;
  delete [] work;
  delete [] iwork;

  return cond;
}

/**
 * API
 *   double* multInvMatMat(CRS TrMat, char tr, double *Matrix)
 *   ��������� ������� �������� � ����������� �� �������
 * INPUT
 *   CRS TrMat - ����������� �������, ������� ���������� 
 *                     �������� � �������� �� �������
 *   char tr         - ��� ����������� �������
 *   double *Matrix  - ���������� ������ �������
 * OUTPUT
 *   
 * RETURN
 *   ���������� ������� 
 **/
double* multInvMatMat(CRS TrMat, char tr, double *Matrix)
{
  double *c;
  char transa = 'N';
  int m = TrMat.N;
  int n = TrMat.N;
  double alpha = 1.0;
  char matdescra[6];
  int ldb = TrMat.N;
  int ldc = TrMat.N;

  matdescra[0] = 'T';
  matdescra[1] = tr;
  matdescra[2] = 'N';
  matdescra[3] = 'C';

  c = new double [TrMat.N * TrMat.N];

  mkl_dcsrsm(&transa, &m, &n, &alpha, matdescra, 
    TrMat.val, TrMat.col_ind, TrMat.row_ptr, TrMat.row_ptr + 1, 
    Matrix, &ldb, c, &ldc);
  return c;
}

/**
 * API
 *   double getConditionNumber(CRS L, CRS U, CRS A)
 *   ������ ����� �������������� � ��������������������
 * INPUT
 *   CRS A - ������� � ������� CRS
 *   CRS L - ������ ������� ������������������� � ������� CRS
 *   CRS U 
 * OUTPUT
 *   
 * RETURN
 *   ����� �������������� 
 **/
double getConditionNumber(CRS L, CRS U, CRS A)
{
    // ��������������� ����������
  double *Matrix, *tmpMatrix;
  int* Size = &(A.N);
  int* ipiv = new int [A.N];
  double* ANorm; 
  double* work = new double [A.N*4];
  int* iwork = new int [A.N];
  int lda = A.N;
  int* Lda = &lda;
  int info = 0;
  int *Info = &info;
  double rcond;
  double* Rcond = &rcond;
  // ����� � ����� ���������������
  double norm, cond = -1.0;

  //�������� ������� � �������� ����
  Matrix = CRStoGeneral(A);

  tmpMatrix = multInvMatMat(L, 'L', Matrix);

  delete [] Matrix;
  Matrix = tmpMatrix;

  tmpMatrix = multInvMatMat(U, 'U', Matrix);

  delete [] Matrix;
  Matrix = tmpMatrix;

  //��������� �����
  norm = MatrixNormI(Matrix, A.N);

  //�LU ����������
  dgetrf(Size,Size,Matrix,Lda,ipiv,Info);
  
  //����� ����� ������ �����
  ANorm = &norm;
  dgecon("O",Size,Matrix,Lda,ANorm,Rcond,work,iwork,Info);
  cond = 1.0/(rcond);

  delete [] Matrix;
  delete [] ipiv;
  delete [] work;
  delete [] iwork;

  return cond;
}