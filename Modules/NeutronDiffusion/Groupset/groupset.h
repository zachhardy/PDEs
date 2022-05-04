#ifndef GROUPSET_H
#define GROUPSET_H

#include "vector.h"
#include "Math/matrix.h"
#include "Math/LinearSolvers/linear_solver.h"

#include <memory>
#include <vector>

namespace neutron_diffusion
{

class Groupset
{
public:
  int id;
  std::vector<size_t> groups;

  /*---------- Options ----------*/
  math::LinearSolverType linear_solver_type = math::LinearSolverType::LU;
  size_t max_iterations = 100;
  double tolerance = 1.0e-8;

  /*---------- System Storage ----------*/
  math::Matrix<double> matrix;
  math::Vector<double> rhs;
  std::shared_ptr<math::LinearSolver<double>> linear_solver;

  Groupset() : id(-1) {}
  explicit Groupset(int groupset_num) : id(groupset_num) {}
};

}

#endif //GROUPSET_H
