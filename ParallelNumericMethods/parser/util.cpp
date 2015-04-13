#include "util.h"

// ������� ���������� ������� � ������� CRS (3 �������, ���������� � ����)
// �������� ������ ��� ���� Value, Col � RowIndex
// ���������� ����� �������� mtx
void InitializeMatrix(int N, int NZ, crsMatrix &mtx)
{
  mtx.N = N;
  mtx.NZ = NZ;
  mtx.Value  = new double[NZ];
  mtx.Col    = new int[NZ];
  mtx.RowIndex = new int[N + 1];
}

// ����������� ������, ���������� ��� ���� mtx
void FreeMatrix(crsMatrix &mtx)
{
  delete [] mtx.Value;
  delete [] mtx.Col;
  delete [] mtx.RowIndex;
  mtx.Value = NULL;
  mtx.Col = NULL;
  mtx.RowIndex = NULL;
}

// ������� ����� imtx � omtx, ������� ������ ��� ���� Value, Col � RowIndex
void CopyMatrix(crsMatrix imtx, crsMatrix &omtx)
{
  // ������������� �������������� �������
  int N  = imtx.N;
  int NZ = imtx.NZ;
  InitializeMatrix(N, NZ, omtx);
  // �����������
  memcpy(omtx.Value   , imtx.Value   , NZ * sizeof(double));
  memcpy(omtx.Col   , imtx.Col   , NZ * sizeof(int));
  memcpy(omtx.RowIndex, imtx.RowIndex, (N + 1) * sizeof(int));
}