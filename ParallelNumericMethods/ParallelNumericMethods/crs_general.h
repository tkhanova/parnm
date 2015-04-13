#pragma once


struct CRS {
	double * val;
	int * col_ind;
	int * row_ptr;

	int N, NZ;
	CRS ( char * filename ) {
		std::ifstream fin ( filename );

		while ( fin.peek() == '%' ) fin.ignore ( 2048, '\n' );

		fin >>  N >> N >> NZ;
		cout << N << " " << NZ << endl;
		val = new double[NZ];
		col_ind = new int[NZ];
		row_ptr = new int[NZ];

		// Read the data in coo
		for ( int i = 0; i < NZ; i++ )  {
			fin >> row_ptr[i] >> col_ind[i] >> val[i];
			cout << row_ptr[i] << " " << col_ind[i] << " " << val[i] << endl;
			row_ptr[i] -= 1;  /* adjust from 1-based to 0-based */
			col_ind[i] -= 1;
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

	CRS() : N ( 0 ), NZ ( 0 ) {}

	~CRS() {
		delete[] val;
		delete[] col_ind;
		delete[] row_ptr;
	}
};

void coo2csr_in ( int n, int nz, double *a, int *i_idx, int *j_idx );
void sort ( int *col_idx, double *a, int start, int end);
int StructTranspose ( int n, int* column, int* row, int* &tColumn, int* &tRow );
void InitializeMatrix(int N, int NZ, CRS &mtx);