#include <algorithm>
#include <cassert>
#include <iostream>
#include <thread>
#include <omp.h>
#include "ProdMatMat.hpp"

namespace {
void prodSubBlocks(int iRowBlkA, int iColBlkB, int iColBlkA, int szBlock,
                   const Matrix& A, const Matrix& B, Matrix& C) {
// Première optimisation : On prend la combinaison la plus optimale
  for (int j = iColBlkB; j < std::min(B.nbCols, iColBlkB + szBlock); j++)
  {
    for (int k = iColBlkA; k < std::min(A.nbCols, iColBlkA + szBlock); k++)
    {
      for (int i = iRowBlkA; i < std::min(A.nbRows, iRowBlkA + szBlock); ++i)
      {
        C(i, j) += A(i, k) * B(k, j);
      }
    }
  }
}
  const int szBlock = 32;
}  // namespace

Matrix operator*(const Matrix& A, const Matrix& B) {
  Matrix C(A.nbRows, B.nbCols, 0.0);
  prodSubBlocks(0, 0, 0, std::max({A.nbRows, B.nbCols, A.nbCols}), A, B, C);
  return C;
}
