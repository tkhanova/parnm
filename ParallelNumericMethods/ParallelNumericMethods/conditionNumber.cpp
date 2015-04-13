#include "conditionNumber.h"
#include "crs_general.h"
#include <math.h>

/**
 * API
 *   double MatrixNormO(double* Matrix, int size)
 *   первая норма матрицы
 * INPUT
 *   double* Matrix - матрица в плотном виде
 *   int     size   - размер матрицы
 * OUTPUT
 *   
 * RETURN
 *   первая норма матрицы
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
 *   бесконечная норма матрицы
 * INPUT
 *   double* Matrix - матрица в плотном виде
 *   int     size   - размер матрицы
 * OUTPUT
 *   
 * RETURN
 *   бесконечная норма матрицы
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
 *   перевод матрицы из формата CRS в плотный
 * INPUT
 *   CRS A - матрица в формате CRS
 * OUTPUT
 *   
 * RETURN
 *   матрица в плотном виде
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
 *   оценка числа обусловлености
 * INPUT
 *   CRS A - матрица в формате CRS
 * OUTPUT
 *   
 * RETURN
 *   число обусловлености 
 **/
double getConditionNumber(CRS A)
{
  // вспомогательные переменные
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
  // нормы и числа обусловленности
  double norm, cond = -1.0;

  //приводим матрицу к плотному виду
  Matrix = CRStoGeneral(A);

  //вычисляем норму
  norm = MatrixNormI(Matrix, A.N);

  //РLU разложение
  dgetrf(Size,Size,Matrix,Lda,ipiv,Info);
  
  //число обусл первая норма
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
 *   умножение матрицы обратной к треугольной на плотную
 * INPUT
 *   CRS TrMat - треугольная матрица, которую необходимо 
 *                     обратить и умножить на плотную
 *   char tr         - тип треугольной матрицы
 *   double *Matrix  - умножаемая справа матрица
 * OUTPUT
 *   
 * RETURN
 *   умноженная матрица 
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
 *   оценка числа обусловлености с предобуславливателем
 * INPUT
 *   CRS A - матрица в формате CRS
 *   CRS L - фактор матрицы предобуславливателя в формате CRS
 *   CRS U 
 * OUTPUT
 *   
 * RETURN
 *   число обусловлености 
 **/
double getConditionNumber(CRS L, CRS U, CRS A)
{
    // вспомогательные переменные
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
  // нормы и числа обусловленности
  double norm, cond = -1.0;

  //приводим матрицу к плотному виду
  Matrix = CRStoGeneral(A);

  tmpMatrix = multInvMatMat(L, 'L', Matrix);

  delete [] Matrix;
  Matrix = tmpMatrix;

  tmpMatrix = multInvMatMat(U, 'U', Matrix);

  delete [] Matrix;
  Matrix = tmpMatrix;

  //вычисляем норму
  norm = MatrixNormI(Matrix, A.N);

  //РLU разложение
  dgetrf(Size,Size,Matrix,Lda,ipiv,Info);
  
  //число обусл первая норма
  ANorm = &norm;
  dgecon("O",Size,Matrix,Lda,ANorm,Rcond,work,iwork,Info);
  cond = 1.0/(rcond);

  delete [] Matrix;
  delete [] ipiv;
  delete [] work;
  delete [] iwork;

  return cond;
}