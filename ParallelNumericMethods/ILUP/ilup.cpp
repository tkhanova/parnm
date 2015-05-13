#include "ilup.h"
#include "queue.h"
#include <omp.h>
#include "validation.h"
using namespace std;



/**
* API
*   int symbolicILUp(int p, int n, int * col, int * row, 
*                    int * lucol, int * lurow, double * &luval,
*                    int * uptr, int &countL, int &countU);
*   символьная фаза ILU(p)
* INPUT
*   int    n     - размер матрицы
*   матрица A
*   int  * col   - индексы колонок матрицы a
*   int  * row   - индексы начала строк матрицы a
*   портрет матрицы LU
*   int  * &lucol - индексы колонок матрицы lu
*   int  * &lurow - индексы начала строк матрицы lu
*   double * &luval - значения матрицы
*   int  * uptr  - индексы диагональных элементов
*                  в массиве luval
* OUTPUT
*   double * luval - значения L и U разложенных матриц
*   int &countL    - размер матрицы L
*   int &countU    - размер матрицы U
* RETURN
*   возвращается код ошибки
**/
int symbolicILUp(int p, int n, int * col, int * row, 
                 int * &lucol, int * &lurow, double * &luval,
                 int * uptr, int &countL, int &countU)
{
  int h, s, f;   // счетчики циклов 
  int jcol;            // и временные переменные
  //queue<int> Q;
  MyQueue Q;
  Q.init(2*p+1); //Магическое число, но вроде больше она не разрастается
  int * len;
  int * adj;
  int * visited;
  len = new int[n];
  adj = new int[n];
  visited = new int[n];
  countL = 0;
  countU = 0;

  for(int j = 0; j < n; j++)
  {
    visited[j] = -1;
    adj[j] = 0;
  }
  int tmp=0;
//#pragma omp parallel for reduction(+:tmp)  private(h, s, jcol, f)
  for(int i = 0; i < n; i++)
  {
	  MyQueue Q;
	  Q.init(2*p+1); 

    Q.push(i);
    len[i] = 0;
    visited[i] = i;
    while(!Q.empty())
    {
      h = Q.pop();
      s = row[h];
      f = row[h + 1];
      for(int j = s; j < f; j++)
      {
        jcol = col[j];
        if(visited[jcol] != i)
        {
//#pragma omp critical
			{visited[jcol] = i;}
          if((jcol > i) && (len[h]<p))
          {
			  Q.push(jcol);
            len[jcol] = len[h] + 1;
          }
          if(jcol < i)
          {
            tmp++;
            adj[i]++;
          }
        }
      }
    }
  }
  countL = tmp;
  for(int i = 0; i < n; i++)
  {
	  s = row[i];
	  f = row[i + 1];
	  int j;
	  for(j = s; (j < f) && (col[j] < i); j++)
		  ;
	  uptr[i] = j;
	  if(col[uptr[i]] != i)
	  {
		  return -(i + 1);
	  }
	  countU += (f - j);
	  adj[i] += (f - j);
  }

  if(lucol != NULL)
  {
	  delete []lucol;
  }
  if(lurow != NULL)
  {
	  delete []lurow;
  }
  if(luval != NULL)
  {
    delete []luval;
  }

  lucol = new int[countL + countU];
  lurow = new int[n + 1];
  luval = new double[countL + countU];

  memset(luval, 0, (countL + countU) * sizeof (double));

  lurow[0] = 0;
  for(int i = 0; i < n; i++)
  {
    lurow[i + 1] = lurow[i] + adj[i];
    //adj[i] = 0;
  }
  memset(adj, 0, n * sizeof (int));

//#pragma omp parallel for  private(h, s, jcol, f)
  for(int i = 0; i < n; i++)
  {
	  MyQueue Q;
	  Q.init(2*p+1); 
	  Q.push(i);
	  len[i] = 0;
//#pragma omp critical
	  {visited[i] = i;}
    while(!Q.empty())
    {
      h = Q.pop();
      s = row[h];
      f = row[h + 1];
      for(int j = s; j < f; j++)
      {
        jcol = col[j];
        if(visited[jcol] != i)
        {
//#pragma omp critical
			{visited[jcol] = i;}
          if((jcol > i) && (len[h]<p))
          {
			  Q.push(jcol);
            len[jcol] = len[h] + 1;
          }
          if(jcol < i)
          {
            lucol[lurow[i] + adj[i]] = jcol;
            adj[i]++;
          }
        }
      }
    }

    s = uptr[i];
    f = row[i + 1];
    uptr[i] = lurow[i] + adj[i];
    for(int j = s; j < f; j++)
    {
      lucol[lurow[i] + adj[i]] = col[j];
      adj[i]++;
    }
  }

  //sort index
  int *tCol;
  int *tRow;
  StructTranspose(n, lucol, lurow, tCol, tRow);
  delete []lucol;
  lucol = NULL;
  delete []lurow;
  lurow = NULL;
  StructTranspose(n, tCol, tRow, lucol, lurow);
  delete []tCol;
  delete []tRow;

  delete[] len;
  delete[] adj;
  delete[] visited;

  return ILU_OK;
}


/**
* API
*   int symbolicILUpWithMultiplication(int p, crsMatrix& A, 
*                    int * lucol, int * lurow, double * &luval,
*                    int * uptr, int &countL, int &countU);
*   символьная фаза ILU(p)
* INPUT
*   int    n     - размер матрицы
*   матрица A
*   int  * col   - индексы колонок матрицы a
*   int  * row   - индексы начала строк матрицы a
*   портрет матрицы LU
*   int  * &lucol - индексы колонок матрицы lu
*   int  * &lurow - индексы начала строк матрицы lu
*   double * &luval - значения матрицы
*   int  * uptr  - индексы диагональных элементов
*                  в массиве luval
* OUTPUT
*   double * luval - значения L и U разложенных матриц
*   int &countL    - размер матрицы L
*   int &countU    - размер матрицы U
* RETURN
*   возвращается код ошибки
**/
int symbolicILUpWithMultiplication(int p, crsMatrix& A, 
	crsMatrix* LU, 
	int * uptr)
{
	//crsMatrix* LU = new crsMatrix;
	crsMatrix* tmp = new crsMatrix;
	int n = A.N;
	if (p==0)
	{
		CopyMatrix(A, *LU);
	}
	else
	{
		CopyMatrix(A, *tmp);
		for (int i=1; i<=p; i++)
		{
			ProductSparseMatrix(A, *tmp, *LU);
			crsMatrix* qqq = tmp;
			tmp = LU;
			LU = qqq;
		}
		delete LU;
		LU = tmp;	
	}

	//fill uptr - index of diagonal elements
	for(int i = 0; i < n; i++)
	{
		int	s = LU->RowIndex[i];
		int f = LU->RowIndex[i + 1];
		int j;
		for(j = s; (j < f) && (LU->Col[j] < i); j++)
			;
		uptr[i] = j;
		
		if(LU->Col[uptr[i]] != i)
		{
			return -(i + 1);
		}
	}

	//lucol = LU->Col;
	//lurow = LU->RowIndex;
//	luval = LU->Value;


	return ILU_OK;
}

/**
* API
*   int numericalILUp(int n, double * a, int * col, int * row, 
*                     int * lucol, int * lurow, int * uptr,
*                     double * luval);
*   численная фаза ILU(p)
* INPUT
*   int    n     - размер матрицы
*   double * a   - не нулевые элементы
*   int  * col   - индексы колонок матрицы a
*   int  * row   - индексы начала строк матрицы a
*   int  * lucol - индексы колонок матрицы lu
*   int  * lurow - индексы начала строк матрицы lu
*   int  * uptr  - индексы диагональных элементов
*                  в массиве luval
* OUTPUT
*   double * luval - значения L и U разложенных матриц
* RETURN
*   возвращается код ошибки
*   0    - разложение выполнено успешно
*   -(n + 1) - номер строки где на диагонале 0
**/
int numericalILUp(int p, int n, double * a, int * col, int * row, 
                  int * lucol, int * lurow, int * uptr,
                  double * luval, int NZ)
{
	bool USE_LEVELS = true;
  int j1, j2;     // граница текущей строки
  int jrow;       // номер текущего столбца
  int k, j, jj;   // счетчики циклов+
  int *iw = NULL;
  int jw;
  double t1;

  iw = new int[n];
  
  int* levels = new int [NZ];

  memset(iw, 0, n * sizeof(int));
  memset(levels, 0, NZ * sizeof(int));

  //копирование исходной матрицы
  j = 0;
  for(k = 0; k < row[n]; k++)
  {
    j1 = col[k];
    j2 = lucol[k + j];
    while(j1 != j2)
    {
      j++;
      j2 = lucol[k + j];
    }
    luval[k + j] = a[k];
  }


  for(k = 0; k < n; k++) //for i = 2 to N do
  {
    j1 = lurow[k];
    j2 = lurow[k + 1];
    for(j = j1; j < j2; j++)
    {
      iw[lucol[j]] = j;
    }
    for(j = j1; (j < j2) && (lucol[j] < k); j++) //for k = 1 to i - 1 and (i, k)  N (| A|p+1) with lev(aik) <= p do
    {
		if (USE_LEVELS)
		{
			if (levels[j]>p) //with lev(aik) <= p
				continue;
		}

      jrow = lucol[j];
      t1 = luval[j] / luval[uptr[jrow]]; //aik = aik/akk
      luval[j] = t1;
	  for(jj = uptr[jrow]+1; jj < lurow[jrow + 1]; jj++) //for j = k + 1 to N and (i, j) ∈ N (| A|p+1)
	  {
        jw = iw[lucol[jj]];
        if(jw != 0)
        {
          luval[jw] = luval[jw] - t1 * luval[jj]; //aij = aij − aik*akj

		  if (USE_LEVELS)
		  {
			  levels[jw] = min(levels[jw], levels[j] + levels[jj] + 1);  //lev(aij ) = min(lev(aij ), lev(aik) + lev(akj ) + 1)
		  }
        }
      }
    }


	if (USE_LEVELS)
		for(j = j1; (j < j2); j++)//for j = 2 to N do
		{
			if (levels[j]>p) //if lev(aij ) > p then
				luval[j] = 0; //aij = 0
		}

    jrow = lucol[j];
    if((jrow != k) || (fabs(luval[j]) < EPSILON))
    {
      break;
    }
    for(j = j1; j < j2; j++)
    {
      iw[lucol[j]] = 0;
    }
  }

  delete [] iw;
  delete[] levels;
  if(k < n)
    return -(k+1);
  return 0;
}