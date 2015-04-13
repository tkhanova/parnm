#include "sparseMatrixOperation.h"

#include <time.h>
#include <memory.h>

/**
 * API
 *   int StructTranspose(int n, int* column, int* row, 
 *                       int* &tColumn, int* &tRow);
 *   ���������������� ��������� �������
 * INPUT
 *   int  n        - ������ �������
 *   int* column   - CRS �������� �������
 *   int* row  
 * OUTPUT
 *   int* &tColumn - CRS �������� �������
 *   int* &tRow  
 * RETURN
 * ������������ ��� ������
**/
int StructTranspose(int n, int* column, int* row, int* &tColumn, int* &tRow)
{
  int i, j;
  int nz;
  int S;
  nz = row[n];

  tColumn = new int [nz];
  tRow    = new int [n + 1];
  
  memset(tRow, 0, (n + 1) * sizeof(int));
  for (i = 0; i < nz; i++) 
  {
    tRow[column[i] + 1]++;
  }

  S = 0;
  for (i = 1; i <= n; i++) 
  {
    int tmp = tRow[i];
    tRow[i] = S;
    S = S + tmp;
  }

  for (i = 0; i < n; i++) 
  {
    int j1 = row[i];
    int j2 = row[i+1];
    int Col = i; // ������� � AT - ������ � �
    for (j = j1; j < j2; j++) 
    {
      int RIndex = column[j];  // ������ � AT
      int IIndex = tRow[RIndex + 1];
      tColumn[IIndex] = Col;
      tRow   [RIndex + 1]++;
    }
  }
  return ILU_OK;
}