#include "petsc_utils.h"

#include "vector.h"
#include "matrix.h"
#include "Sparse/sparse_matrix.h"


using namespace Math;


void
PETScUtils::CreateSquareMatrix(Mat& A, PetscInt n)
{
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n);
  MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
  MatSetFromOptions(A);
  MatSetUp(A);
}


void
PETScUtils::CreateMatrix(Mat& A, const Matrix& mat)
{
  CreateSquareMatrix(A, static_cast<PetscInt>(mat.n_rows()));

  for (PetscInt i = 0; i < mat.n_rows(); ++i)
    for (PetscInt j = 0; j < mat.n_cols(); ++j)
    {
      PetscScalar a_ij = static_cast<PetscScalar>(mat(i, j));
      if (a_ij!= 0.0)
        MatSetValue(A, i, j, a_ij, INSERT_VALUES);
    }

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
}



void
PETScUtils::CreateMatrix(Mat& A, const SparseMatrix& mat)
{
  CreateSquareMatrix(A, static_cast<PetscInt>(mat.n_rows()));

  // Add data to the PETSc matrix
  for (const auto el : mat)
    MatSetValue(A, el.row, el.column, el.value, ADD_VALUES);

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
}
