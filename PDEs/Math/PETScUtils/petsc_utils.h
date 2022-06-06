#ifndef PETSC_UTILS_H
#define PETSC_UTILS_H

#include <petscksp.h>
#include <vector>


namespace Math
{
  //########## Forward declarations
  class Vector;
  class Matrix;
  class SparseMatrix;
}


namespace Math::PETScUtils
{
  void CreateVector(Vec& x, PetscInt n);
  void CreateVector(Vec& x, const Vector& vec);
  void CreateVector(Vec& x, const std::vector<double>& vec);

  void CreateSquareMatrix(Mat& A, PetscInt n);
  void CreateMatrix(Mat& A, const Matrix& mat);
  void CreateMatrix(Mat& A, const SparseMatrix& mat);

//  void InitMatrixSparsity(Mat& A, PetscInt n,
//                          const std::vector<PetscInt>& nnz_per_row);

  void CopyToVector(Vec& x, Vector& vec);
  void CopyToSTLVector(Vec& x, std::vector<double>& vec);
}

#endif //PETSC_UTILS_H
