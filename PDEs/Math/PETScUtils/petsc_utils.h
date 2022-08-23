#ifndef PETSC_UTILS_H
#define PETSC_UTILS_H

#include <petscksp.h>
#include <vector>


namespace PDEs
{
  namespace Math
  {
    // forward declarations
    class Vector;
    class Matrix;
    class SparseMatrix;


    namespace PETScUtils
    {
      /**
       * Create a PETSc vector of size \p n
       */
      void
      CreateVector(Vec& x, PetscInt n);

      /**
       * Create a PETSc vector from a Vector.
       */
      void
      CreateVector(Vec& x, const Vector& vec);

      /**
       * Create a PETSc vector from and STL vector of doubles.
       */
      void
      CreateVector(Vec& x, const std::vector<double>& vec);

      /**
       * Create a square PETSc matrix of dimension \p n.
       */
      void
      CreateMatrix(Mat& A, PetscInt n);

      /**
       * Create a PETSc matrix with \p n_rows and \p n_cols.
       */
      void
      CreateMatrix(Mat& A, PetscInt n_rows, PetscInt n_cols);

      /**
       * Create a PETSc matrix from a Matrix.
       */
      void
      CreateMatrix(Mat& A, const Matrix& mat);

      /**
       * Create a PETSc matrix from a SparseMatrix.
       */
      void
      CreateMatrix(Mat& A, const SparseMatrix& mat);

      //  void InitMatrixSparsity(Mat& A, PetscInt n,
      //                          const std::vector<PetscInt>& nnz_per_row);

      /**
       * Copy the contents of a PETSc vector \p x to Vector \p vec.

       \param[in] x A PETSc vector to copy from.
       \param[out] vec A Vector to copy to.
       */
      void
      CopyToVector(Vec& x, Vector& vec);

      /**
       * Copy the contents of a PETSc vector \p x to an STL vector \p vec.
       *
       * \param[in] x A PETSc vector to copy from.
       * \param[out] vec An STL vector to copy to.
       */
      void
      CopyToVector(Vec& x, std::vector<double>& vec);
    }
  }
}


#endif //PETSC_UTILS_H
