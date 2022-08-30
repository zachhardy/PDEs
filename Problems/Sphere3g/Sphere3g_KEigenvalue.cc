#include "ortho_grids.h"

#include "material.h"
#include "CrossSections/cross_sections.h"

#include "timer.h"

#include "Math/LinearSolvers/Iterative/cg.h"
#include "Math/LinearSolvers/Direct/cholesky.h"
#include "LinearSolvers/PETSc/petsc_solver.h"

#include "NeutronDiffusion/KEigenvalueSolver/keigenvalue_solver.h"

#include <iostream>
#include <vector>


using namespace PDEs;
using namespace Grid;
using namespace Math;
using namespace Physics;
using namespace LinearSolvers;
using namespace NeutronDiffusion;


int main(int argc, char** argv)
{
  //============================================================
  // Mesh
  //============================================================

  size_t n_cells = 100;
  double slab_width = 6.0;
  double cell_width = slab_width / (double) n_cells;

  std::vector<double> vertices(1, 0.0);
  for (size_t i = 0; i < n_cells; ++i)
    vertices.emplace_back(vertices.back() + cell_width);

  auto coord_sys = CoordinateSystemType::SPHERICAL;
  auto mesh = create_1d_orthomesh(vertices, coord_sys);

  //============================================================
  // Materials
  //============================================================

  auto material = std::make_shared<Material>();

  // Create the cross sections
  auto xs = std::make_shared<CrossSections>();
  xs->read_xs_file("Problems/Sphere3g/xs/base3g.xs", 0.05);
  material->properties.emplace_back(xs);

  // Create the multigroup source
  std::vector<double> mg_source(xs->n_groups, 1.0);
  auto src = std::make_shared<IsotropicMultiGroupSource>(mg_source);
  material->properties.emplace_back(src);

  //============================================================
  // Linear Solver
  //============================================================

  Options opts;
  opts.verbosity = 0;
  opts.tolerance = 1.0e-10;
  opts.max_iterations = 10000;

  std::shared_ptr<LinearSolverBase<SparseMatrix>> linear_solver;
  linear_solver = std::make_shared<PETScSolver>(KSPCG, PCLU, opts);

  //============================================================
  // Create the diffusion solver
  //============================================================

  KEigenvalueSolver solver;

  solver.mesh = mesh;
  solver.materials.emplace_back(material);
  solver.linear_solver = linear_solver;

  solver.verbosity = 1;
  solver.use_precursors = false;

  solver.algorithm = Algorithm::DIRECT;

  //============================================================
  // Define boundary conditions
  //============================================================

  solver.boundary_info.emplace_back(BoundaryType::REFLECTIVE, -1);
  solver.boundary_info.emplace_back(BoundaryType::ZERO_FLUX, -1);

  //============================================================
  // Run the problem
  //============================================================

  PetscInitialize(&argc, &argv, (char*) 0, NULL);

  Timer timer;

  solver.initialize();

  timer.start();
  solver.execute();
  timer.stop();

  PetscFinalize();

  std::cout << timer.get_time() << " ms\n";
}
