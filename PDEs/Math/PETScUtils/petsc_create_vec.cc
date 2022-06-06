#include "petsc_utils.h"
#include "vector.h"


using namespace Math;


void
PETScUtils::CreateVector(Vec& x, PetscInt n)
{
  VecCreate(PETSC_COMM_WORLD, &x);
  VecSetSizes(x, PETSC_DECIDE, n);
  VecSetFromOptions(x);
  VecSetUp(x);
}


void
PETScUtils::CreateVector(Vec& x, const Vector& vec)
{
  CreateVector(x, vec.size());

  const double* v_ptr = vec.data();
  for (size_t i = 0; i < vec.size(); ++i)
    VecSetValue(x, i, *v_ptr++, INSERT_VALUES);
}


void
PETScUtils::CreateVector(Vec& x, const std::vector<double>& vec)
{
  CreateVector(x, vec.size());

  const double* v_ptr = vec.data();
  for (size_t i = 0; i < vec.size(); ++i)
    VecSetValue(x, i, *v_ptr++, INSERT_VALUES);
}
