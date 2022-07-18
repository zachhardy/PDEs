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

  size_t n_x = 21, n_y = 21;
  double X = 80.0, Y = 80.0;
  double dx = X / (n_x - 1), dy = Y / (n_y - 1);

  std::vector<double> x_verts(1, 0.0);
  for (size_t i = 0; i < n_x - 1; ++i)
    x_verts.push_back(x_verts.back() + dx);

  std::vector<double> y_verts(1, 0.0);
  for (size_t i = 0; i < n_y - 1; ++i)
    y_verts.push_back(y_verts.back() + dy);

  auto mesh = create_2d_orthomesh(x_verts, y_verts, true);

  for (auto& cell : mesh->cells)
  {
    auto& c = cell.centroid;
    if (24.0 < c.x < 56.0 and 24.0 < c.y < 56.0)
      cell.material_id = 0;
    else if (0.0 < c.x < 24.0 and 24.0 < c.y < 56.0)
      cell.material_id = 1;
    else if (24.0 < c.x < 56.0 and 0.0 < c.y < 24.0)
      cell.material_id = 1;
    else
      cell.material_id = 2;
  }

  //============================================================
  // Materials
  //============================================================
  using namespace Physics;

  std::vector<std::shared_ptr<Material>> materials;
  materials.emplace_back(std::make_shared<Material>("Fuel 0"));
  materials.emplace_back(std::make_shared<Material>("Fuel 1"));
  materials.emplace_back(std::make_shared<Material>("Fuel 2"));

  std::vector<std::shared_ptr<CrossSections>> xs;
  for (size_t i = 0; i < materials.size(); ++i)
    xs.emplace_back(std::make_shared<CrossSections>());

  std::vector<std::string> xs_paths;
  xs_paths.emplace_back("Test/TWIGL/xs/mat0.xs");
  xs_paths.emplace_back("Test/TWIGL/xs/mat0.xs");
  xs_paths.emplace_back("Test/TWIGL/xs/mat1.xs");

  for (size_t i = 0; i < materials.size(); ++i)
  {
    xs[i]->read_xs_file(xs_paths[i]);
    materials[i]->properties.emplace_back(xs[i]);
  }

  const size_t n_groups = xs.front()->n_groups;

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

  for (auto& material : materials)
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
  solver.boundary_info.emplace_back(BoundaryType::ZERO_FLUX, -1);
  solver.boundary_info.emplace_back(BoundaryType::ZERO_FLUX, -1);
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
