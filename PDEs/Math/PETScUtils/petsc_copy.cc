#include "petsc_utils.h"

#include "vector.h"
#include <vector>


using namespace Math;


void
PETScUtils::CopyToVector(Vec& x, Vector& vec)
{
  PetscInt size;
  VecGetSize(x, &size);

  vec.clear();
  vec.resize(size);

  // Get pointer to data to copy to
  double* v_ptr = vec.data();

  // Get pointer to PETSc data
  const double* x_ptr;
  VecGetArrayRead(x, &x_ptr);

  // Perform the copying
  const double* end_ptr = x_ptr + size;
  while(x_ptr != end_ptr)
    *v_ptr++ = *x_ptr++;

  // Remove access to PETSc data
  VecRestoreArrayRead(x, &x_ptr);
}


void
PETScUtils::CopyToVector(Vec& x, std::vector<double>& vec)
{
  PetscInt size;
  VecGetSize(x, &size);

  vec.clear();
  vec.resize(size);

  // Get pointer to data to copy to
  double* v_ptr = vec.data();

  // Get pointer to PETSc data
  const double* x_ptr;
  VecGetArrayRead(x, &x_ptr);

  // Perform the copying
  const double* end_ptr = x_ptr + size;
  while(x_ptr != end_ptr)
    *v_ptr++ = *x_ptr++;

  // Remove access to PETSc data
  VecRestoreArrayRead(x, &x_ptr);
}
