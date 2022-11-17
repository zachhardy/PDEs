#include "ortho_grids.h"

#include "material.h"
#include "CrossSections/cross_sections.h"

#include "timer.h"

#include "Math/LinearSolvers/Iterative/cg.h"
#include "Math/LinearSolvers/Direct/cholesky.h"
#include "LinearSolvers/PETSc/petsc_solver.h"

#include "NeutronDiffusion/SteadyStateSolver/steadystate_solver.h"
#include "NeutronDiffusion/KEigenvalueSolver/keigenvalue_solver.h"
#include "NeutronDiffusion/TransientSolver/transient_solver.h"

#include <iomanip>
#include <iostream>
#include <vector>


using namespace std;
using namespace PDEs;
using namespace Grid;
using namespace Math;
using namespace Physics;
using namespace LinearSolvers;
using namespace NeutronDiffusion;


int main(int argc, char** argv)
{
  double radius = 6.0;
  double density = 0.05;
  double sigs_01 = 1.46;

  string xsdir = "xs";
  string outdir = "outputs";

  for (int i = 0; i < argc; ++i)
  {
    string arg(argv[i]);
    cout << "Parsing argument " << i << " " << arg << endl;

    if (arg.find("radius") == 0)
      radius = stod(arg.substr(arg.find('=') + 1));
    else if (arg.find("density") == 0)
      density = stod(arg.substr(arg.find('=') + 1));
    else if (arg.find("scatter") == 0)
      sigs_01 = stod(arg.substr(arg.find('=') + 1));
    else if (arg.find("output_directory") == 0)
      outdir = arg.substr(arg.find('=') + 1);
    else if (arg.find("xs_directory") == 0)
      xsdir = arg.substr(arg.find('=') + 1);
  }

  //============================================================
  // Mesh
  //============================================================

  size_t n_cells = 100;
  double cell_width = radius / (double) n_cells;

  vector<double> vertices(1, 0.0);
  for (size_t i = 0; i < n_cells; ++i)
    vertices.emplace_back(vertices.back() + cell_width);

  auto coord_sys = CoordinateSystemType::SPHERICAL;
  auto mesh = create_1d_orthomesh(vertices, coord_sys);

  //============================================================
  // Materials
  //============================================================

  auto material = make_shared<Material>();

  // Create the cross sections
  auto xs = make_shared<CrossSections>();
  xs->read_xs_file(xsdir+"/base3g.xs", density);
  xs->transfer_matrices[0][1][0] = sigs_01 * density;
  material->properties.emplace_back(xs);

  //============================================================
  // Linear Solver
  //============================================================

  Options opts;
  opts.verbosity = 0;
  opts.tolerance = 1.0e-10;
  opts.max_iterations = 10000;

  shared_ptr<LinearSolverBase<SparseMatrix>> linear_solver;
  linear_solver = make_shared<PETScSolver>(KSPCG, PCLU, opts);

  //============================================================
  // Create the diffusion solver
  //============================================================

  TransientSolver solver;

  solver.mesh = mesh;
  solver.materials.emplace_back(material);
  solver.linear_solver = linear_solver;

  solver.verbosity = 0;
  solver.use_precursors = false;

  solver.algorithm = Algorithm::DIRECT;

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
      { return 1.0 - p.z() * p.z() / (radius * radius); };
  solver.initial_conditions[0] = ic;
  solver.initial_conditions[1] = ic;

  solver.adaptive_time_stepping = false;
  solver.coarsen_threshold = 0.01;
  solver.refine_threshold = 0.025;

  //============================================================
  // Run the problem
  //============================================================

  PetscInitialize(&argc, &argv, (char*) 0, NULL);

  Timer timer;

  timer.start();
  solver.initialize();
  solver.execute();
  timer.stop();

  PetscFinalize();

  cout << "\nSimulation Time: " << setprecision(10)
            << timer.get_time() << " ms\n";
}
