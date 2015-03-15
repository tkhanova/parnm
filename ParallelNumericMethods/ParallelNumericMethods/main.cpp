#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <queue>

using namespace std;

void sort ( int *col_idx, double *a, int start, int end ) {
    int i, j, it;
    double dt;
    
    for ( i = end - 1; i > start; i-- )
        for ( j = start; j < i; j++ )
            if ( col_idx[j] > col_idx[j + 1] ) {
                if ( a ) {
                    dt = a[j];
                    a[j] = a[j + 1];
                    a[j + 1] = dt;
                }
                
                it = col_idx[j];
                col_idx[j] = col_idx[j + 1];
                col_idx[j + 1] = it;
            }
}


/* converts COO format to CSR format, in-place,
if SORT_IN_ROW is defined, each row is sorted in column index.
On return, i_idx contains row_start position */

void coo2csr_in ( int n, int nz, double *a, int *i_idx, int *j_idx ) {
    int *row_start;
    int i, j;
    int init, i_next, j_next, i_pos;
    double dt, a_next;
    row_start = ( int * ) malloc ( ( n + 1 ) * sizeof ( int ) );
    
    if ( !row_start ) {
        printf ( "coo2csr_in: cannot allocate temporary memory\n" );
        exit ( 1 );
    }
    
    for ( i = 0; i <= n; i++ ) row_start[i] = 0;
    
    /* determine row lengths */
    for ( i = 0; i < nz; i++ ) row_start[i_idx[i] + 1]++;
    
    for ( i = 0; i < n; i++ ) row_start[i + 1] += row_start[i];
    
    for ( init = 0; init < nz; ) {
        dt = a[init];
        i = i_idx[init];
        j = j_idx[init];
        i_idx[init] = -1;
        
        while ( 1 ) {
            i_pos = row_start[i];
            a_next = a[i_pos];
            i_next = i_idx[i_pos];
            j_next = j_idx[i_pos];
            a[i_pos] = dt;
            j_idx[i_pos] = j;
            i_idx[i_pos] = -1;
            row_start[i]++;
            
            if ( i_next < 0 ) break;
            
            dt = a_next;
            i = i_next;
            j = j_next;
        }
        
        init++;
        
        while ( ( i_idx[init] < 0 ) && ( init < nz ) )  init++;
    }
    
    /* shift back row_start */
    for ( i = 0; i < n; i++ ) i_idx[i + 1] = row_start[i];
    
    i_idx[0] = 0;
    
    for ( i = 0; i < n; i++ ) {
        sort ( j_idx, a, i_idx[i], i_idx[i + 1] );
    }
}


struct CRS {
    double * val;
    int * col_ind;
    int * row_ptr;
    
    int N, M, NZ;
    
    CRS ( char * filename ) {
        std::ifstream fin ( filename );
        
        while ( fin.peek() == '%' ) fin.ignore ( 2048, '\n' );
        
        fin >> M >> N >> NZ;
        cout << M << " " << N << " " << NZ << endl;
        val = new double[NZ];
        col_ind = new int[NZ];
        row_ptr = new int[N + 1];
        
        // Read the data in coo
        for ( int i = 0; i < NZ; i++ )  {
            fin >> row_ptr[i] >> col_ind[i] >> val[i];
            row_ptr[i] -= 1;  /* adjust from 1-based to 0-based */
            col_ind[i] -= 1;
            cout << row_ptr[i] << " " << col_ind[i] << " " << val[i] << endl;
        }
        
        fin.close();
        coo2csr_in ( N, NZ, val, row_ptr, col_ind );
        
        for ( int i = 0; i < NZ; i++ )  {
            cout << val[i] << " ";
        }
        
        cout << endl;
        
        for ( int i = 0; i < NZ; i++ )  {
            cout << col_ind[i] << " ";
        }
        
        cout << endl;
        
        for ( int i = 0; i < N + 1; i++ )  {
            cout << row_ptr[i] << " ";
        }
        
        cout << endl;
    }
    
    CRS() {}
};

/** * API * int symbolicILUp(int p, int n, int * col, int * row, * int * lucol, int * lurow, double * &luval, * int * uptr, int &countL, int &countU);
* символическая фаза ILU(p) * INPUT * int n - размер матрицы * матрица A * int * col - индексы колонок матрицы a * int * row - индексы начала строк матрицы a
* портрет матрицы LU
* int * &lucol - индексы колонок матрицы lu
* int * &lurow - индексы начала строк матрицы lu * double * &luval - значения матрицы
* int * uptr - индексы диагональных элементов * в массиве luval
* OUTPUT
* double * luval - значения L и U разложенных матриц
* int &countL - размер матрицы L
* int &countU - размер матрицы U
* RETURN * возвращается код ошибки **/

int symbolicILUp ( int p, int n, int * col, int * row, int * &lucol, int * &lurow, double * &luval, int * uptr ) {
    int i, j, h, s, f; // счетчики циклов
    int jcol; // и временные переменные
    int * len;
    int * visited;
    len = new int[n];
    int * adj;
    adj = new int[n];
    visited = new int[n];
    int countL = 0;
    int countU = 0;
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
    //int *tCol;
    //int *tRow;
    //StructTranspose(n, lucol, lurow, tCol, tRow);
    //delete []lucol;
    //delete []lurow;
    //StructTranspose(n, tCol, tRow, lucol, lurow);
    //delete []tCol;
    //delete []tRow;
    delete[] len;
    delete[] adj;
    return 0;
}

int numericalILUp ( int n, double * a, int * col, int * row, int * lucol, int * lurow, int * uptr, double * luval ) {
    int j1, j2; // граница текущей строки
    int jrow; // номер текущего столбца
    int k, j, jj; // счетчики циклов
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

void ILUp ( CRS & matrix, CRS & LU ) {
    int * upptr = new int[matrix.N];
    symbolicILUp ( matrix.NZ, matrix.N, matrix.col_ind, matrix.row_ptr, LU.col_ind, LU.row_ptr, LU.val, upptr );
    numericalILUp ( matrix.N, matrix.val, matrix.col_ind, matrix.row_ptr, LU.col_ind,  LU.row_ptr, upptr, LU.val );
}


int main ( int argv, char *argc[] ) {
    if ( argv < 2 ) return 1;
    
    CRS matrix ( argc[1] );
    CRS LU;
    ILUp ( matrix, LU );
    return 0;
}