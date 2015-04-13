#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <queue>
#include <mkl.h>
#include "conditionNumber.h"
#include "crs_general.h"

using namespace std;


/** * API * int symbolicILUp(int p, int n, int * col, int * row, * int * lucol, int * lurow, double * &luval, * int * uptr, int &countL, int &countU);
* ������������� ���� ILU(p) * INPUT * int n - ������ ������� * ������� A * int * col - ������� ������� ������� a * int * row - ������� ������ ����� ������� a
* ������� ������� LU
* int * &lucol - ������� ������� ������� lu
* int * &lurow - ������� ������ ����� ������� lu * double * &luval - �������� �������
* int * uptr - ������� ������������ ��������� * � ������� luval
* OUTPUT
* double * luval - �������� L � U ����������� ������
* int &countL - ������ ������� L
* int &countU - ������ ������� U
* RETURN * ������������ ��� ������ **/

int symbolicILUp ( int p, int n, int * col, int * row, int * &lucol, int * &lurow, double * &luval, int * uptr, int & countL, int & countU ) {
    int i, j, h, s, f; // �������� ������
    int jcol; // � ��������� ����������
    int * len;
    int * visited;
    len = new int[n];
    int * adj;
    adj = new int[n];
    visited = new int[n];
    countL = 0;
    countU = 0;
    queue<int> Q;
    
    for ( j = 0; j < n; j++ ) {
        visited[j] = -1;
        adj[j] = 0;
    }
    
    for ( i = 0; i < n; i++ ) {
        Q.push ( i );
        len[i] = 0;
        visited[i] = i;
        
        while ( ! ( Q.empty() ) ) {
            h = Q.front();
            Q.pop();
            s = row[h];
            f = row[h + 1];
            
            for ( j = s; j < f; j++ ) {
                jcol = col[j];
                
                if ( visited[jcol] != i ) {
                    visited[jcol] = i;
                    
                    if ( ( jcol > i ) && ( len[h] < p ) ) {
                        Q.push ( jcol );
                        len[jcol] = len[h] + 1;
                    }
                    
                    if ( jcol < i ) {
                        countL++;
                        adj[i]++;
                    }
                }
            }
        }
    }
    
    for ( i = 0; i < n; i++ ) {
        s = row[i];
        f = row[i + 1];
        
        for ( j = s; ( j < f ) && ( col[j] < i ); j++ );
        
        uptr[i] = j;
        
        if ( col[uptr[i]] != i ) {
            return - ( i + 1 );
        }
        
        countU += ( f - j );
        adj[i] += ( f - j );
    }
    
    //if ( lucol != NULL ) {
    //    delete []lucol;
    //}
    //
    //if ( lurow != NULL ) {
    //    delete []lurow;
    //}
    //
    //if ( luval != NULL ) {
    //    delete []luval;
    //}
    lucol = new int[countL + countU];
    lurow = new int[n + 1];
    luval = new double[countL + countU];
    memset ( luval, 0, ( countL + countU ) * sizeof ( double ) );
    lurow[0] = 0;
    
    for ( i = 0; i < n; i++ ) {
        lurow[i + 1] = lurow[i] + adj[i];
        adj[i] = 0;
    }
    
    for ( i = 0; i < n; i++ ) {
        Q.push ( i );
        len[i] = 0;
        visited[i] = i;
        
        while ( ! ( Q.empty() ) ) {
            h = Q.front();
            Q.pop();
            s = row[h];
            f = row[h + 1];
            
            for ( j = s; j < f; j++ ) {
                jcol = col[j];
                
                if ( visited[jcol] != i ) {
                    visited[jcol] = i;
                    
                    if ( ( jcol > i ) && ( len[h] < p ) ) {
                        Q.push ( jcol );
                        len[jcol] = len[h] + 1;
                    }
                    
                    if ( jcol < i ) {
                        lucol[lurow[i] + adj[i]] = jcol;
                        adj[i]++;
                    }
                }
            }
        }
        
        s = uptr[i];
        f = row[i + 1];
        uptr[i] = lurow[i] + adj[i];
        
        for ( j = s; j < f; j++ ) {
            lucol[lurow[i] + adj[i]] = col[j];
            adj[i]++;
        }
    }
    
    //sort index
    int *tCol;
    int *tRow;
    StructTranspose ( n, lucol, lurow, tCol, tRow );
    delete []lucol;
    delete []lurow;
    StructTranspose ( n, tCol, tRow, lucol, lurow );
    delete []tCol;
    delete []tRow;
    delete[] len;
    delete[] adj;
    return 0;
}

int numericalILUp ( int n, double * a, int * col, int * row, int * lucol, int * lurow, int * uptr, double * luval ) {
    int j1, j2; // ������� ������� ������
    int jrow; // ����� �������� �������
    int k, j, jj; // �������� ������
    int *iw = NULL;
    int jw;
    double t1;
    iw = new int[n];
    memset ( iw, 0, n * sizeof ( int ) );
    j = 0;
    
    for ( k = 0; k < row[n]; k++ ) {
        j1 = col[k];
        j2 = lucol[k + j];
        
        while ( j1 != j2 ) {
            j++;
            j2 = lucol[k + j];
        }
        
        luval[k + j] = a[k];
    }
    
    for ( k = 0; k < n; k++ ) {
        j1 = lurow[k];
        j2 = lurow[k + 1];
        
        for ( j = j1; j < j2; j++ ) {
            iw[lucol[j]] = j;
        }
        
        for ( j = j1; ( j < j2 ) && ( lucol[j] < k ); j++ ) {
            jrow = lucol[j];
            t1 = luval[j] / luval[uptr[jrow]];
            luval[j] = t1;
            
            for ( jj = uptr[jrow] + 1; jj < lurow[jrow + 1]; jj++ ) {
                jw = iw[lucol[jj]];
                
                if ( jw != 0 ) {
                    luval[jw] = luval[jw] - t1 * luval[jj];
                }
            }
        }
        
        jrow = lucol[j];
        
        if ( ( jrow != k ) || ( fabs ( luval[j] ) < 1e-8 ) ) {
            break;
        }
        
        for ( j = j1; j < j2; j++ ) {
            iw[lucol[j]] = 0;
        }
    }
    
    delete [] iw;
    
    if ( k < n ) return - ( k + 1 );
    
    return 0;
}


/**
 * API
 *  void ProductSparseMatrix (CRS &A,  CRS &B, 
 *                            CRS &C);
 *  ��������� ����������� ������ C = A * B
 * INPUT
 *   CRS &A   - ������� A
 *   CRS &B   - ������� B
 * OUTPUT
 *   CRS &C   - C = A * B
 * RETURN
 */
void ProductSparseMatrix(CRS &A, CRS &B, CRS &C)
{
  int n = A.N;

  // �������� ��������� ��� ������ ������� MKL
  // ��������������� ������� A � B � �������
  int i, j;
  for (i = 0; i < A.NZ; i++)
    A.col_ind[i]++;
  for (i = 0; i < B.NZ; i++)
    B.col_ind[i]++;
  for (j = 0; j <= n; j++)
  {
    A.row_ptr[j]++;
    B.row_ptr[j]++;
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
  C.row_ptr = new int[n + 1];
// ��������� ���������� ��������� ��������� � ������� C
  request = 1;
  C.val = 0;
  C.col_ind = 0;
  mkl_dcsrmultcsr(&trans, &request, &sort, &n, &n, &n, A.val, A.col_ind, 
                  A.row_ptr, B.val, B.col_ind, B.row_ptr, C.val, C.col_ind,
                  C.row_ptr, &nzmax, &info);
  int nzc = C.row_ptr[n] - 1;
  C.val = new double[nzc];
  C.col_ind   = new int[nzc];
// ��������� C = A * B
  request = 2;
  mkl_dcsrmultcsr(&trans, &request, &sort, &n, &n, &n, A.val, A.col_ind,
                  A.row_ptr, B.val, B.col_ind, B.row_ptr, C.val, C.col_ind,
                  C.row_ptr, &nzmax, &info);
  C.N = n;
  C.NZ = nzc;

  // �������� � ����������� ���� ������� A, B � �
  for (i = 0; i < A.NZ; i++)
    A.col_ind[i]--;
  for (i = 0; i < B.NZ; i++)
    B.col_ind[i]--;
  for (i = 0; i < C.NZ; i++)
    C.col_ind[i]--;
  for (j = 0; j <= n; j++)
  {
    A.row_ptr[j]--;
    B.row_ptr[j]--;
    C.row_ptr[j]--;
  }
}

/**
 * API
 *  double MatrixCompare (CRS &A,  CRS &M);
 *  ������� ������� ������� ������
 * INPUT
 *   CRS &A   - �������� �������
 *   CRS &M   - �������������������
 * OUTPUT
 *   
 * RETURN
 *  ������� �������
 */
double MatrixCompare (CRS &A,  CRS &M)
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
    i  = M.row_ptr[k];
    fM = M.row_ptr[k + 1];
    j  = A.row_ptr[k];
    fA = A.row_ptr[k + 1];
    while((i < fM) && (j<fA))
    {
      if(M.col_ind[i] != A.col_ind[j])
      {
        i++;
      } 
      else 
      {
        if(diff < fabs(M.val[i] - A.val[j]))
        {
          diff =  fabs(M.val[i] - A.val[j]);
        }
        j++;
        i++;
      }
    }
  }
  return diff; 
}

/**
 * API
 *  void LUmatrixSeparation (CRS ilu, int *uptr, 
 *                           CRS &L, CRS &U);
 *  ���������� ������� �� L � U �������
 * INPUT
 *   CRS ilu  - ������� ilu � ����� ���������
 *   int    * uptr  - ������� ������������ ���������
 *                    � ������� �������� ilu
 * OUTPUT
 *   CRS &L   - ���������� ������� L
 *   CRS &U   - ���������� ������� U
 * RETURN
 */
void LUmatrixSeparation(CRS ilu, int *uptr, CRS &L, CRS &U)
{
  int countL, countU;
  int i, j, s, f, k;
  double *val;
  int    *col;
  countU = 0;
  for(i = 0; i < ilu.N; i++)
  {
    countU += (ilu.row_ptr[i+1] - uptr[i]);
  }
  countL = ilu.NZ + ilu.N - countU;
  InitializeMatrix(ilu.N, countL, L);
  InitializeMatrix(ilu.N, countU, U);
  // fill L matrix
  k = 0;
  val = L.val;
  col = L.col_ind;
  L.row_ptr[0] = k;
  for(i = 0; i < ilu.N; i++)
  {
    s = ilu.row_ptr[i];
    f = uptr[i];
    for(j = s; j < f; j++)
    {
      val[k] = ilu.val[j];
      col[k] = ilu.col_ind[j];
      k++;
    }
    val[k] = 1.0;
    col[k] = i;
    k++;
    L.row_ptr[i + 1] = k;
  }
  // fill U matrix
  k = 0;
  val = U.val;
  col = U.col_ind;
  U.row_ptr[0] = k;
  for(i = 0; i < ilu.N; i++)
  {
    s = uptr[i];
    f = ilu.row_ptr[i + 1];
    for(j = s; j < f; j++)
    {
      val[k] = ilu.val[j];
      col[k] = ilu.col_ind[j];
      k++;
    }
    U.row_ptr[i + 1] = k;
  }
}


void ILUp ( int p, CRS & matrix, CRS & LU , int*& upptr) {
    //int * upptr = new int[matrix.N];
    int countL, countU;
    symbolicILUp ( p, matrix.N, matrix.col_ind, matrix.row_ptr, LU.col_ind, LU.row_ptr, LU.val, upptr, countL, countU );
    LU.N = matrix.N;
    LU.NZ = countL + countU;
    numericalILUp ( matrix.N, matrix.val, matrix.col_ind, matrix.row_ptr, LU.col_ind,  LU.row_ptr, upptr, LU.val );
}


int main ( int argv, char *argc[] ) {
    if ( argv < 2 ) return 1;
    
    CRS matrix ( argc[1] );
    cout << matrix << endl;
    CRS LU;
    int p = 2;  // should be less than N
	int* upptr = new int[matrix.N];
    ILUp ( p, matrix, LU, upptr);

	CRS L,U,M;

	LUmatrixSeparation(LU, upptr, L, U);

	ProductSparseMatrix(L, U, M);
	// ���������� ������� ������� ��������
	double diff = MatrixCompare(matrix, M);

	printf("distinction value of matrix %f \n", diff);

	if(matrix.N < 4000)
	{
		double condA, condMA;
		condA = getConditionNumber(*matA);
		printf("Condition number of matrix        A: %lf\n", condA);
		condMA = getConditionNumber(L, U, *matA);
		printf("Condition number of matrix M^-1 * A: %lf\n", condMA);
	}

    cout << LU << endl;
    return 0;
}