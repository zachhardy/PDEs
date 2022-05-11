#ifndef GROUPSET_H
#define GROUPSET_H

#include "vector.h"
#include "matrix.h"
#include "LinearSolvers/linear_solver.h"

#include <memory>
#include <vector>
#include <cinttypes>


namespace neutron_diffusion
{

class Groupset
{
public:
  int id;
  std::vector<uint64_t> groups;

  /*---------- Options ----------*/
  math::LinearSolverType linear_solver_type = math::LinearSolverType::LU;
  uint64_t max_iterations = 100;
  double tolerance = 1.0e-8;

  /*---------- System Storage ----------*/
  math::SparseMatrix<double> matrix;
  math::Vector<double> rhs;

  Groupset() : id(-1) {}
  explicit Groupset(int groupset_num) : id(groupset_num) {}
};

}

#endif //GROUPSET_H
