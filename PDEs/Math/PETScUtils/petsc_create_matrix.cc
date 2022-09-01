#include "petsc_utils.h"

#include "vector.h"
#include "matrix.h"
#include "Math/sparse_matrix.h"


using namespace PDEs;
using namespace Math;


void
PETScUtils::create_petsc_matrix(Mat& A, PetscInt n)
{ create_petsc_matrix(A, n, n); }


void
PETScUtils::create_petsc_matrix(Mat& A, PetscInt n_rows, PetscInt n_cols)
{
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n_rows, n_cols);
  MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
  MatSetFromOptions(A);
  MatSetUp(A);
}


void
PETScUtils::create_petsc_matrix(Mat& A, const Matrix& mat)
{
  create_petsc_matrix(A, static_cast<PetscInt>(mat.n_rows()),
                      static_cast<PetscInt>(mat.n_cols()));

  for (PetscInt i = 0; i < mat.n_rows(); ++i)
    for (PetscInt j = 0; j < mat.n_cols(); ++j)
    {
      PetscScalar a_ij = static_cast<PetscScalar>(mat(i, j));
      if (a_ij != 0.0)
        MatSetValue(A, i, j, a_ij, INSERT_VALUES);
    }

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
}


void
PETScUtils::create_petsc_matrix(Mat& A, const SparseMatrix& mat)
{
  create_petsc_matrix(A, static_cast<PetscInt>(mat.n_rows()),
                      static_cast<PetscInt>(mat.n_cols()));

  // Add data to the PETSc matrix
  for (const auto el: mat)
    MatSetValue(A, el.row, el.column, el.value, ADD_VALUES);

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
}
