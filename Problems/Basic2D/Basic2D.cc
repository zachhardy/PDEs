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
  //============================================================
  // Mesh
  //============================================================
  using namespace Grid;

  size_t n_x = 3, n_y = 3;
  double X = 1.0, Y = 1.0;
  double dx = X / (n_x - 1), dy = Y / (n_y - 1);

  std::vector<double> x_verts(1, 0.0);
  for (size_t i = 0; i < n_x - 1; ++i)
    x_verts.push_back(x_verts.back() + dx);

  std::vector<double> y_verts(1, 0.0);
  for (size_t i = 0; i < n_y - 1; ++i)
    y_verts.push_back(y_verts.back() + dy);

  auto mesh = create_2d_orthomesh(x_verts, y_verts, true);

  //============================================================
  // Materials
  //============================================================
  using namespace Physics;

  auto material = std::make_shared<Material>();

  auto xs = std::make_shared<CrossSections>();
  xs->read_xs_file("Problems/Basic2D/xs/test_1g.xs");
  material->properties.emplace_back(xs);

  const size_t n_groups = xs->n_groups;

  //============================================================
  // Linear Solver
  //============================================================
  using namespace LinearSolver;

  Options opts;
  opts.verbosity = 0;
  opts.tolerance = 1.0e-10;
  opts.max_iterations = 10000;

  std::shared_ptr<LinearSolverBase<SparseMatrix>> linear_solver;
  linear_solver = std::make_shared<PETScSolver>(KSPCG, PCLU, opts);

  //============================================================
  // Create the diffusion solver
  //============================================================
  using namespace NeutronDiffusion;
  KEigenvalueSolver solver;
  solver.mesh = mesh;
  solver.materials.emplace_back(material);
  solver.linear_solver = linear_solver;

  solver.verbosity = 1;
  solver.use_precursors = true;

  solver.solution_technique = SolutionTechnique::FULL_SYSTEM;

  //============================================================
  // Initialize groups and groupsets
  //============================================================

  for (size_t g = 0; g < n_groups; ++g)
    solver.groups.emplace_back(g);

  Groupset groupset(0);
  for (size_t g = 0; g < n_groups; ++g)
    groupset.groups.emplace_back(solver.groups[g]);
  solver.groupsets.emplace_back(groupset);

  //============================================================
  // Define boundary conditions
  //============================================================

  solver.boundary_info.emplace_back(BoundaryType::REFLECTIVE, -1);
  solver.boundary_info.emplace_back(BoundaryType::REFLECTIVE, -1);
  solver.boundary_info.emplace_back(BoundaryType::REFLECTIVE, -1);
  solver.boundary_info.emplace_back(BoundaryType::REFLECTIVE, -1);

  //============================================================
  // Run the problem
  //============================================================

  PetscInitialize(&argc,&argv,(char*)0,NULL);

  Timer timer;

  solver.initialize();

  timer.start();
  solver.execute();
  timer.stop();

  PetscFinalize();

  std::cout << timer.get_time() << " ms\n";
}
