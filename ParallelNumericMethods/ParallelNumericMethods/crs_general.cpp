#include "crs_general.h"
#include <memory>

using namespace std;


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

/**
 * API
 *   int StructTranspose(int n, int* column, int* row,
 *                       int* &tColumn, int* &tRow);
 *   Транспонирование структуры матрицы
 * INPUT
 *   int  n        - размер матрицы
 *   int* column   - CRS описание матрицы
 *   int* row
 * OUTPUT
 *   int* &tColumn - CRS описание матрицы
 *   int* &tRow
 * RETURN
 * возвращается код ошибки
**/
int StructTranspose ( int n, int* column, int* row, int* &tColumn, int* &tRow ) {
    int i, j;
    int nz;
    int S;
    nz = row[n];
    tColumn = new int [nz];
    tRow    = new int [n + 1];
    memset ( tRow, 0, ( n + 1 ) * sizeof ( int ) );
    
    for ( i = 0; i < nz; i++ ) {
        tRow[column[i] + 1]++;
    }
    
    S = 0;
    
    for ( i = 1; i <= n; i++ ) {
        int tmp = tRow[i];
        tRow[i] = S;
        S = S + tmp;
    }
    
    for ( i = 0; i < n; i++ ) {
        int j1 = row[i];
        int j2 = row[i + 1];
        int Col = i; // Столбец в AT - строка в А
        
        for ( j = j1; j < j2; j++ ) {
            int RIndex = column[j];  // Строка в AT
            int IIndex = tRow[RIndex + 1];
            tColumn[IIndex] = Col;
            tRow   [RIndex + 1]++;
        }
    }
    
    return 0;
}




// Создает квадратную матрицу в формате CRS (3 массива, индексация с нуля)
// Выделяет память под поля Value, Col и row_ptr
// Возвращает через параметр mtx
void InitializeMatrix(int N, int NZ, CRS &mtx)
{
	mtx.N = N;
	mtx.NZ = NZ;
	mtx.val  = new double[NZ];
	mtx.col_ind    = new int[NZ];
	mtx.row_ptr = new int[N + 1];
}
ostream & operator<< ( ostream & s, const CRS & m ) {
    for ( int i = 0; i < m.NZ; i++ )
        s << m.val[i] << " ";
        
    s << endl;
    
    for ( int i = 0; i < m.NZ; i++ )
        s << m.col_ind[i] << " ";
        
    s << endl;
    
    for ( int i = 0; i < m.N + 1; i++ )
        s << m.row_ptr[i] << " ";
        
    s << endl;
    return s;
}