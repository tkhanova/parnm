#include "readMTX.h"

/**
* API
*   int ReadMTXFile(FILE* matrixFile, int* n, int* countNotZero, int** columns, 
*          int** rows, FLOAT_TYPE** values, int** countValueInRow,
*          int* isOrdered, int* isStoreOnlyUpperTriangle, 
*          int* isStoreOnlyLowTriangle)
*   ������ ������� �� mtx ������� ����� � ������������ ������ �������� ������
* INPUT
*   FILE* matrixFile   - ���� �������
* OUTPUT
*   int  n         - ������ �������
*   int* countNotZero  - ���-�� ��������� ���������
*   int** columns    - ������� � ������������ �������
*   int** rows
*   FLOAT_TYPE** values
*   int** countValueInRow - ���������� �������� � ������
*   int* isOrdered    - ����������� �� �������
*   int* isStoreOnlyUpperTriangle - �������� �� ������ ������� �����������
*   int* isStoreOnlyLowTriangle   - �������� �� ������ ������ �����������
* RETURN
*   ������������ ��� ������
**/
int ReadMTXFile(FILE* matrixFile, int* n, int* countNotZero, int** columns, 
        int** rows, FLOAT_TYPE** values, int** countValueInRow,
        int* isOrdered, int* isStoreOnlyUpperTriangle, 
        int* isStoreOnlyLowTriangle)
{
  int error = ILU_OK;
  int i;
  //������ ��� ������ ������ �� �����
  char* sBuf0 = 0;
  char* sBuf1 = 0;
  char* sBuf2 = 0;
  char* sBuf3 = 0;
  int sizeStringBuffer = 256;
  int countUseBuffer;
  int countColumn, countRow;
  size_t charSize = sizeof(char);
  size_t intSize = sizeof(int);

  int currerntColumn;
  int currentRow;
  FLOAT_TYPE currentValue;

  int isFirstValueZero = 1;
  int* minColumnInRow = 0;


  sBuf0 = (char*)malloc(sizeStringBuffer * charSize);
  sBuf1 = (char*)malloc(sizeStringBuffer * charSize);
  sBuf2 = (char*)malloc(sizeStringBuffer * charSize);
  sBuf3 = (char*)malloc(sizeStringBuffer * charSize);
  for (i = 0; i < sizeStringBuffer; i++)
  {
  sBuf0[i] = 0;
  sBuf1[i] = 0;
  sBuf2[i] = 0;
  sBuf3[i] = 0;
  }
  (*n) = -1;

  while ((sBuf0 != 0) && ((*n) <= 0))
  {
  sBuf0 = fgets(sBuf0, sizeStringBuffer, matrixFile);
  countUseBuffer = sscanf(sBuf0, "%s %s %s", sBuf1, sBuf2, sBuf3);
  countRow = atoi(sBuf1);
  countColumn = atoi(sBuf2);
  (*countNotZero) = atoi(sBuf3);
  if ((countRow > 0) && (countColumn > 0))
  {
    if (countRow !=  countColumn)
    {
    error = MATRIX_NOT_SQUARTE;
    break;
    }
    else
    {
    (*n) = countColumn;
    }
  }
  }
  if (sBuf0 == 0)
  {
  error = NO_MATRIX_IN_FILE;
  }
  if (error != 0)
  {
  return error;
  }

  (*countValueInRow) = (int*)malloc((*n) * intSize);
  minColumnInRow = (int*)malloc((*n) * intSize);


  (*columns) = (int*)malloc((*countNotZero) * intSize);
  (*rows) = (int*)malloc((*countNotZero) * intSize);
  (*values) = (FLOAT_TYPE*)malloc((*countNotZero) * sizeof(FLOAT_TYPE));

  for (i = 0; i < (*n); i++)
  {
  (*countValueInRow)[i] = 0;
  minColumnInRow[i] = 0;
  }

  //������ �����
  for (i = 0; i < (*countNotZero); i++)
  {
  countUseBuffer = fscanf(matrixFile, "%s %s %s", sBuf0, sBuf1 , sBuf2);
  if (countUseBuffer == 3)
  {
    currentRow = atoi(sBuf0);
    currerntColumn = atoi(sBuf1);
    currentValue = atof(sBuf2);

    currentRow--;
    currerntColumn--;

    (*rows)[i] = currentRow;
    (*columns)[i] = currerntColumn;
    (*values)[i] = currentValue;

    (*countValueInRow)[currentRow]++;
    if (minColumnInRow[currentRow] >= currerntColumn)
    {
    minColumnInRow[currentRow] = currerntColumn;
    if ((*countValueInRow)[currentRow] > 1)
    {
      (*isOrdered) = 0;
    }
    }
    if (currerntColumn < currentRow)
    {
    (*isStoreOnlyUpperTriangle) = 0;
    }

    if (currerntColumn > currentRow)
    {
    (*isStoreOnlyLowTriangle) = 0;
    }
    if ((currentRow == 0) && (currerntColumn == 0))
    {
    isFirstValueZero = 0;
    }

  }
  else
  {
    error = MISMATCH_MATRIX_FORMAT;
  }
  }

  if (isFirstValueZero == 1)
  {
  error = FIRST_VALUE_ZERO;
  }


  if (sBuf0 != 0)
  {
  free(sBuf0);
  }
  if (sBuf1 != 0)
  {
  free(sBuf1);
  }
  if (sBuf2 != 0)
  {
  free(sBuf2);
  }
  if (sBuf3 != 0)
  {
  free(sBuf3);
  }
  if (minColumnInRow != 0)
  {
  free(minColumnInRow);
  }

  return error;
}

/**
* API
*   CheckSymmetric(int countNotZero, int* columns, int* rows, 
*          FLOAT_TYPE* values)
*   �������� ������� �� ��������������
* INPUT
*   int countNotZero   - ����� ��������� ��������� � ������� 
* OUTPUT
*   int* columns     - ������� � ������������ �������
*   int* rows
*   FLOAT_TYPE* values
* RETURN
*   ������������ ����������� ������� ��� ���
**/
int CheckSymmetric(int countNotZero, int* columns, int* rows, 
           FLOAT_TYPE* values)
{
  int i, j;
  int isHaveSymmetricValue = 1;
  int isSymmetricMatrix = 1;

  for(i = 0; i < countNotZero; i++)
  {
  isHaveSymmetricValue = 1;
  for (j = 0; j < countNotZero; j++)
  {
    if ((fabs(values[i] - values[j]) > EPSILON) && 
    (columns[i] == rows[j]) && 
    (rows[i] == columns[j]) && 
    (i != j))
    {
    isHaveSymmetricValue = 0;
    break;
    }
  }
  if(!isHaveSymmetricValue)
  {
    isSymmetricMatrix = 0;
    break;
  }
  }

  return isSymmetricMatrix;
}


/**
* API
*   SortMatrix(int n, int* column, int* row, FLOAT_TYPE* val);
*   �������������� �������
* INPUT
*   int  n         - ������ �������
*   int* column      - CRS �������� �������
*   int* row  
*   FLOAT_TYPE* val
* OUTPUT
*   int* column      - CRS �������� �������
*   int* row  
*   FLOAT_TYPE* val
* RETURN
* ������������ ����� ������
**/
void SortMatrix(int n, int* column, int* row, FLOAT_TYPE* val)
{
  int i, j;
  int nz;
  int s;
  int j1; 
  int j2; 
  int col;
  double V;
  int RIndex;
  int IIndex;
  int tmp;
  int* tColumn; 
  int* tRow; 
  FLOAT_TYPE* tVal;
  nz = row[n];
  s = 0;
//AllocMatrix(n, nz, &tColumn, &tRow, &tVal);
  tColumn = new int[ nz ];
  tRow = new int[ n + 1 ];
  tVal = new FLOAT_TYPE[ nz ];


  memset(tRow, 0, (n + 1) * sizeof(int));
  for (i = 0; i < nz; i++) 
  tRow[column[i] + 1]++;

  for (i = 1; i <= n; i++) 
  {
  tmp = tRow[i];
  tRow[i] = s;
  s = s + tmp;
  }

  for (i = 0; i < n; i++) 
  {
  j1 = row[i];
  j2 = row[i+1];
  col = i; // ������� � AT - ������ � �
  for (j = j1; j < j2; j++) 
  {
    V = val[j];  // ��������
    RIndex = column[j];  // ������ � AT
    IIndex = tRow[RIndex + 1];
    tVal[IIndex] = V;
    tColumn  [IIndex] = col;
    tRow[RIndex + 1]++;
  }
  }

  memset(row, 0, (n + 1) * sizeof(int));
  s = 0;
  for (i = 0; i < nz; i++) 
  row[tColumn[i] + 1]++;

  for (i = 1; i <= n; i++) 
  {
  tmp = row[i];
  row[i] = s;
  s = s + tmp;
  }

  for (i = 0; i < n; i++) 
  {
  j1 = tRow[i];
  j2 = tRow[i+1];
  col = i; // ������� � AT - ������ � �
  for (j = j1; j < j2; j++) 
  {
    V = tVal[j];  // ��������
    RIndex = tColumn[j];  // ������ � AT
    IIndex = row[RIndex + 1];
    val[IIndex] = V;
    column  [IIndex] = col;
    row[RIndex + 1]++;
  }
  }  

  //FreeMatrixMTX(&tColumn, &tRow, &tVal);
  delete[] tColumn;
  delete[] tRow;
  delete[] tVal;
}

/**
* API
*   ReadMatrixFromFile(char* matrixName, int* n, int** column, 
*            int** row, FLOAT_TYPE** val)
*   ������ ������� �� mtx ������� ����� � crs ������ �������� ������
* INPUT
*   char* matrixName
* OUTPUT
*   int  n         - ������ �������
*   int** column     - CRS �������� �������
*   int** row
*   FLOAT_TYPE** val
* RETURN
*   ������������ ��� ������
**/
int ReadMatrixFromFile(char* matrixName, int* n, int* nz, int** column, 
             int** row, FLOAT_TYPE** val, int *typeOfMatrix)
{
  int i;
  int countNotZero = -1;//����� �� ������� ��������� � ������
  int countNotZeroOld;
  size_t intSize = sizeof(int);
  FILE* matrixFile  = 0;   
  int error = ILU_OK;
  int* countValueInRow = 0;
  int* currentIndexInRow = 0;
  int* columns = 0;
  int* rows = 0;
  FLOAT_TYPE* values = 0; 
  int isSymmetricMatrix = 1;
  int isStoreOnlyUpperTriangle = 1;
  int isStoreOnlyLowTriangle = 1;
  int isOrdered = 1;
  int transposeBuffer = 0;
  int countNotZeroInUpperTriangle;
  int countNotZeroInLowTriangle;
  int countNotZeroInDiagonal;
  int isDiagonalZero = 0;

  fopen_s(&matrixFile, matrixName, "r");
  if(matrixFile == 0)
  {
  printf("Error open files %s\n", matrixName);
  return CANT_OPEN_FILE;
  }

  (*typeOfMatrix) = UNDEFINE_TYPE;

  //������ ����
  error = ReadMTXFile(matrixFile, n, &countNotZero, &columns, &rows, &values, 
  &countValueInRow, &isOrdered, &isStoreOnlyUpperTriangle, 
  &isStoreOnlyLowTriangle);

  if (error != 0)
  {
  fclose(matrixFile);
  return error;
  }

  countNotZeroOld = countNotZero;
  currentIndexInRow = (int*)malloc((*n) * intSize);
  for (i = 0; i < (*n); i++)
  {
  currentIndexInRow[i] = 0;
  }

  if(isStoreOnlyUpperTriangle == 1)
  {
  (*typeOfMatrix) = UPPER_TRIANGULAR;
  }
  //���� ������� ������ ������ �������������, �� �������������
  if (isStoreOnlyLowTriangle == 1)
  {
  isStoreOnlyLowTriangle = 0;
  isStoreOnlyUpperTriangle = 1;
  isOrdered = 0;
  for (i = 0; i < (*n); i++)
  {
    countValueInRow[i] = 0;
  }
  for(i = 0; i < countNotZero; i++)
  {
    transposeBuffer = rows[i];
    rows[i] = columns[i];
    columns[i] = transposeBuffer;
    countValueInRow[rows[i]]++;
  }
  (*typeOfMatrix) = UPPER_TRIANGULAR;
  }


  //���� ������� ������ �� ������ ������� ����������� � ���������, 
  //�� ��������� �� ��������������
  if (!isStoreOnlyUpperTriangle)
  {
  isSymmetricMatrix = CheckSymmetric(countNotZero, columns, rows, values);

  countNotZeroInUpperTriangle = 0;
  countNotZeroInLowTriangle = 0;
  countNotZeroInDiagonal = 0;    

  //���� ������� ������������ �������� ������ ����������� �������
  if (isSymmetricMatrix == 1)
  {
    (*typeOfMatrix) = SYMMETRIC;
  }
  }

  //��������� ������ ��� �������
  //AllocMatrix((*n), countNotZero, column, row, val);
  *column = new int[ countNotZero ];
  *row = new int[ (*n) + 1 ];
  *val = new FLOAT_TYPE[ countNotZero ];
  
  for (i = 0; i < countNotZero; i++)
  {
  (*column)[i] = -1;
  (*val)[i] = -1;
  }
  
  for (i = 0; i <= (*n); i++)
  {
  (*row)[i] = -1;
  }

  (*row)[0] = 0;
  //���������� �������� ������ ����� �������
  for (i = 0; i < (*n); i++)
  {
  (*row)[i + 1] = (*row)[i] + countValueInRow[i];
  }

  //���������� ������ ������� � �������� �������
  for (i = 0; i < countNotZeroOld; i++)
  {
  if (rows[i] != -1)
  {
    (*column)[(*row)[rows[i]] + currentIndexInRow[rows[i]]] = columns[i];
    (*val)[(*row)[rows[i]] + currentIndexInRow[rows[i]]] = values[i];
    currentIndexInRow[rows[i]]++;    
  }
  }

  if (isDiagonalZero == 1)
  {
  for (i = 0; i < (*n); i++)
  {
    if((*column)[(*row)[i]] == -1)
    {
    (*column)[(*row)[i]] = i;
    (*val)[(*row)[i]] = 0;
    }
  }
  }

  //���� ������� �� �������������, �� �������������
  if (!isOrdered)
  {
  SortMatrix(*n, *column, *row, *val);
  }
  if (countValueInRow != 0)
  {
  free(countValueInRow);
  }
  if (currentIndexInRow != 0)
  {
  free(currentIndexInRow);
  }
  if (columns != 0)
  {
  free(columns);
  }
  if (rows != 0)
  {
  free(rows);
  }
  if (values != 0)
  {
  free(values);
  }
  if (matrixFile != 0)
  {
  fclose(matrixFile);
  }
  *nz = countNotZero;
  return error;
}