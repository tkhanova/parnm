#include "routines.h"

void ParseArgv( int argc, char** argv, char* &matrixFile )
{
	matrixFile = argv[ 1 ];
}

void PrintMatrixInCRSView( const crsMatrix* const A )
{
	printf( "Value = " );
	for( int i = 0; i < A->NZ; i++ )
	{
		printf( "%lf ", A->Value[ i ] );
	}
	printf( "\n" );

	printf( "Col = " );
	for( int i = 0; i < A->NZ; i++ )
	{
		printf( "%d ", A->Col[ i ] );
	}
	printf( "\n" );
	
	printf( "RowIndex = " );
	for( int i = 0; i < (A->N + 1); i++ )
	{
		printf( "%d ", A->RowIndex[ i ] );
	}
	printf( "\n\n" );

}

void PrintMatrixInNormalView( const crsMatrix* const A )
{
	int position = 0;
	for( int i = 0; i < A->N; i++ )
	{
		for( int j = 0; j < A->N; j++ )
		{
			if( (A->Col[ position ] == j) && (position < A->RowIndex[ i + 1]) )
			{
				printf( "%10.2f ", A->Value[ position ] );
				position++;
			}
			else
			{
				printf( "%10.2lf ", 0.0 );
			}
		}
		printf( "\n" );
	}
	printf( "\n" );
}


crsMatrix* UpTriangleMatrixToFullSymmetricMatrix( const crsMatrix* const A )
{
	crsMatrix* FullA;
	int countElementsInRowA = 0;
	int* countElementsInRows = new int[ A->N ];
	int* currentPostitionsInRow = new int[ A->N ];

	FullA = new crsMatrix;

	InitializeMatrix( A->N, (A->NZ * 2 - A->N), *FullA);
	
	for( int i = 0; i < A->N; i++ )
	{
		countElementsInRows[ i ] = 0;
		currentPostitionsInRow[ i ] = 0;
	}

	for( int i = 0; i < A->N; i++ )
	{
		for( int k = A->RowIndex[ i ] + 1; k < A->RowIndex[ i + 1 ]; k++ )
		{
			countElementsInRows[ A->Col[ k ] ]++;
		}
	}
	FullA->RowIndex[ 0 ] = 0;
	for( int i = 1; i <= A->N; i++ )
	{
		countElementsInRowA = A->RowIndex[ i ] - A->RowIndex[ i - 1 ];
		FullA->RowIndex[ i ] = FullA->RowIndex[ i - 1 ] + countElementsInRowA + countElementsInRows[ i - 1 ];
		currentPostitionsInRow[ i - 1 ] = FullA->RowIndex[ i - 1 ];
	}
	
	for( int i = 0; i < A->N; i++ )
	{
		for( int j = A->RowIndex[ i ] + 1; j < A->RowIndex[ i + 1 ] ; j++ )
		{
			FullA->Value[ currentPostitionsInRow[ A->Col[ j ] ] ] = A->Value[ j ];
			FullA->Col[ currentPostitionsInRow[ A->Col[ j ] ] ] = i;
			currentPostitionsInRow[ A->Col[ j ] ]++;
		}
	}

	for( int i = 0; i < A->N; i++ )
	{
		for( int j = A->RowIndex[ i ]; j < A->RowIndex[ i + 1 ]; j++ )
		{
			FullA->Value[ currentPostitionsInRow[ i ] ] = A->Value[ j ];
			FullA->Col[ currentPostitionsInRow[ i ] ] = A->Col[ j ];
			currentPostitionsInRow[ i ]++;
		}
	}

	delete[] currentPostitionsInRow;
	delete[] countElementsInRows;
	return FullA;
}

