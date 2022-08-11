#include "ortho_grids.h"

#include "material.h"
#include "CrossSections/cross_sections.h"

#include "timer.h"

#include "LinearSolvers/DirectSolvers"
#include "LinearSolvers/IterativeSolvers"
#include "LinearSolvers/PETSc/petsc_solver.h"

#include "NeutronDiffusion/groupset.h"
#include "NeutronDiffusion/SteadyStateSolver/steadystate_solver.h"
#include "NeutronDiffusion/KEigenvalueSolver/keigenvalue_solver.h"
#include "NeutronDiffusion/TransientSolver/transient_solver.h"

#include <iostream>
#include <vector>


int main(int argc, char** argv)
{
  double radius = 6.0;
  double density = 0.05;
  double sigs_01 = 1.46;
  std::string outdir = "Problems/Sphere3g/outputs";

  for (int i = 0; i < argc; ++i)
  {
    std::string arg(argv[i]);
    std::cout << "Parsing argument " << i << " " << arg << std::endl;

    if (arg.find("radius") == 0)
      radius = std::stod(arg.substr(arg.find("=") + 1));
    else if (arg.find("density") == 0)
      density = std::stod(arg.substr(arg.find("=") + 1));
    else if (arg.find("scatter") == 0)
      sigs_01 = std::stod(arg.substr(arg.find("=") + 1));
    else if (arg.find("output_directory") == 0)
      outdir = arg.substr(arg.find("=") + 1);
  }

  //============================================================
  // Mesh
  //============================================================
  using namespace Grid;

  size_t n_cells = 100;
  double cell_width = radius/(double)n_cells;

  std::vector<double> vertices(1, 0.0);
  for (size_t i = 0; i < n_cells; ++i)
    vertices.emplace_back(vertices.back() + cell_width);

  auto coord_sys = CoordinateSystemType::SPHERICAL;
  auto mesh = create_1d_orthomesh(vertices, coord_sys);

  //============================================================
  // Materials
  //============================================================
  using namespace Physics;

  auto material = std::make_shared<Material>();

  // Create the cross sections
  auto xs = std::make_shared<CrossSections>();
  xs->read_xs_file("Problems/Sphere3g/xs/base3g.xs", density);
  xs->transfer_matrices[0][1][0] = sigs_01 * density;
  material->properties.emplace_back(xs);

  size_t n_groups = xs->n_groups;

  // Create the multigroup source
  std::vector<double> mg_source(n_groups, 1.0);
  auto src = std::make_shared<IsotropicMultiGroupSource>(mg_source);
  material->properties.emplace_back(src);

  //============================================================
  // Linear Solver
  //============================================================
  using namespace Math;
  using namespace Math::LinearSolver;

  Options opts;
  opts.verbosity = 0;
  opts.tolerance = 1.0e-10;
  opts.max_iterations = 10000;

  std::shared_ptr<LinearSolverBase<SparseMatrix>> linear_solver;
//    linear_solver = std::make_shared<CG>(opts);
  linear_solver = std::make_shared<PETScSolver>(KSPCG, PCLU, opts);

  //============================================================
  // Create the diffusion solver
  //============================================================
  using namespace NeutronDiffusion;

  TransientSolver solver;
  solver.mesh = mesh;
  solver.materials.emplace_back(material);
  solver.linear_solver = linear_solver;

  solver.verbosity = 0;
  solver.use_precursors = false;

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

  //============================================================
  // Define transient parameters
  //============================================================

  solver.t_end = 0.1;
  solver.dt = solver.t_end / 50;
  solver.time_stepping_method = TimeSteppingMethod::CRANK_NICHOLSON;
  solver.normalization_method = NormalizationMethod::TOTAL_POWER;

  solver.write_outputs = true;
  solver.output_directory = outdir;

  auto ic = [radius](const Point p)
  { return 1.0 - p.z*p.z/(radius*radius); };
  solver.initial_conditions[0] = ic;
  solver.initial_conditions[1] = ic;

  solver.adaptivity = false;
  solver.coarsen_threshold = 0.01;
  solver.refine_threshold = 0.025;

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
