#include "ilup.h"

using namespace std;

/**
* API
*   int symbolicILUp(int p, int n, int * col, int * row, 
*                    int * lucol, int * lurow, double * &luval,
*                    int * uptr, int &countL, int &countU);
*   ���������� ���� ILU(p)
* INPUT
*   int    n     - ������ �������
*   ������� A
*   int  * col   - ������� ������� ������� a
*   int  * row   - ������� ������ ����� ������� a
*   ������� ������� LU
*   int  * &lucol - ������� ������� ������� lu
*   int  * &lurow - ������� ������ ����� ������� lu
*   double * &luval - �������� �������
*   int  * uptr  - ������� ������������ ���������
*                  � ������� luval
* OUTPUT
*   double * luval - �������� L � U ����������� ������
*   int &countL    - ������ ������� L
*   int &countU    - ������ ������� U
* RETURN
*   ������������ ��� ������
**/
int symbolicILUp(int p, int n, int * col, int * row, 
                 int * &lucol, int * &lurow, double * &luval,
                 int * uptr, int &countL, int &countU)
{
  int i, j, h, s, f;   // �������� ������ 
  int jcol;            // � ��������� ����������
  queue<int> Q;
  int * len;
  int * adj;
  int * visited;
  len = new int[n];
  adj = new int[n];
  visited = new int[n];
  countL = 0;
  countU = 0;

  for(j = 0; j < n; j++)
  {
    visited[j] = -1;
    adj[j] = 0;
  }
  for(i = 0; i < n; i++)
  {
    Q.push(i);
    len[i] = 0;
    visited[i] = i;
    while(!(Q.empty()))
    {
      h = Q.front();
      Q.pop();
      s = row[h];
      f = row[h + 1];
      for(j = s; j < f; j++)
      {
        jcol = col[j];
        if(visited[jcol] != i)
        {
          visited[jcol] = i;
          if((jcol > i) && (len[h]<p))
          {
            Q.push(jcol);
            len[jcol] = len[h] + 1;
          }
          if(jcol < i)
          {
            countL++;
            adj[i]++;
          }
        }
      }
    }
  }

  for(i = 0; i < n; i++)
  {
    s = row[i];
    f = row[i + 1];
    for(j = s; (j < f) && (col[j] < i); j++);
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
  for(i = 0; i < n; i++)
  {
    lurow[i + 1] = lurow[i] + adj[i];
    adj[i] = 0;
  }

  for(i = 0; i < n; i++)
  {
    Q.push(i);
    len[i] = 0;
    visited[i] = i;
    while(!(Q.empty()))
    {
      h = Q.front();
      Q.pop();
      s = row[h];
      f = row[h + 1];
      for(j = s; j < f; j++)
      {
        jcol = col[j];
        if(visited[jcol] != i)
        {
          visited[jcol] = i;
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
    for(j = s; j < f; j++)
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
  delete []lurow;
  StructTranspose(n, tCol, tRow, lucol, lurow);
  delete []tCol;
  delete []tRow;

  delete[] len;
  delete[] adj;

  return ILU_OK;
}

/**
* API
*   int numericalILUp(int n, double * a, int * col, int * row, 
*                     int * lucol, int * lurow, int * uptr,
*                     double * luval);
*   ��������� ���� ILU(p)
* INPUT
*   int    n     - ������ �������
*   double * a   - �� ������� ��������
*   int  * col   - ������� ������� ������� a
*   int  * row   - ������� ������ ����� ������� a
*   int  * lucol - ������� ������� ������� lu
*   int  * lurow - ������� ������ ����� ������� lu
*   int  * uptr  - ������� ������������ ���������
*                  � ������� luval
* OUTPUT
*   double * luval - �������� L � U ����������� ������
* RETURN
*   ������������ ��� ������
*   0    - ���������� ��������� �������
*   -(n + 1) - ����� ������ ��� �� ��������� 0
**/
int numericalILUp(int n, double * a, int * col, int * row, 
                  int * lucol, int * lurow, int * uptr,
                  double * luval)
{
  int j1, j2;     // ������� ������� ������
  int jrow;       // ����� �������� �������
  int k, j, jj;   // �������� ������+
  int *iw = NULL;
  int jw;
  double t1;

  iw = new int[n];

  memset(iw, 0, n * sizeof(int));

  //����������� �������� �������
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


  for(k = 0; k < n; k++)
  {
    j1 = lurow[k];
    j2 = lurow[k + 1];
    for(j = j1; j < j2; j++)
    {
      iw[lucol[j]] = j;
    }
    for(j = j1; (j < j2) && (lucol[j] < k); j++)
    {
      jrow = lucol[j];
      t1 = luval[j] / luval[uptr[jrow]];
      luval[j] = t1;
      for(jj = uptr[jrow]+1; jj < lurow[jrow + 1]; jj++)
      {
        jw = iw[lucol[jj]];
        if(jw != 0)
        {
          luval[jw] = luval[jw] - t1 * luval[jj];
        }
      }
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
  if(k < n)
    return -(k+1);
  return 0;
}