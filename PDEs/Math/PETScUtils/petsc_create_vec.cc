#include "petsc_utils.h"
#include "vector.h"


using namespace PDEs;
using namespace Math;


void
PETScUtils::create_petsc_vector(Vec& x, PetscInt n)
{
  VecCreate(PETSC_COMM_WORLD, &x);
  VecSetSizes(x, PETSC_DECIDE, n);
  VecSetFromOptions(x);
  VecSetUp(x);
}


void
PETScUtils::create_petsc_vector(Vec& x, const Vector& vec)
{
  create_petsc_vector(x, vec.size());

  const double* v_ptr = vec.data();
  for (size_t i = 0; i < vec.size(); ++i)
    VecSetValue(x, i, *v_ptr++, INSERT_VALUES);
}


void
PETScUtils::create_petsc_vector(Vec& x, const std::vector<double>& vec)
{
  create_petsc_vector(x, vec.size());

  const double* v_ptr = vec.data();
  for (size_t i = 0; i < vec.size(); ++i)
    VecSetValue(x, i, *v_ptr++, INSERT_VALUES);
}
