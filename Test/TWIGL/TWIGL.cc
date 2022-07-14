#include "mesh.h"
#include "ortho_grids.h"

#include "material.h"
#include "CrossSections/cross_sections.h"

#include "timer.h"

#include "LinearSolvers/DirectSolvers"
#include "LinearSolvers/IterativeSolvers"
#include "LinearSolvers/PETSc/petsc_solver.h"

#include "NeutronDiffusion/groupset.h"
#include "NeutronDiffusion/KEigenvalueSolver/keigenvalue_solver.h"

#include <iostream>
#include <vector>
#include <map>

#include <petsc.h>


int main(int argc, char** argv)
{
  using namespace Grid;

  size_t n_x = 41, n_y = 41;
  double X = 80.0, Y = 80.0;
  double dx = X / (n_x - 1), dy = Y / (n_y - 1);

  std::vector<double> x_verts(1, 0.0);
  for (size_t i = 0; i < n_x - 1; ++i)
    x_verts.push_back(x_verts.back() + dx);

  std::vector<double> y_verts(1, 0.0);
  for (size_t i = 0; i < n_y - 1; ++i)
    y_verts.push_back(y_verts.back() + dy);

  auto mesh = create_2d_orthomesh(x_verts, y_verts, true);
}
