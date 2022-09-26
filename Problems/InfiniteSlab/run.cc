#include "ortho_grids.h"

#include "material.h"
#include "CrossSections/cross_sections.h"

#include "timer.h"

#include "Math/LinearSolvers/Iterative/cg.h"
#include "Math/LinearSolvers/Direct/cholesky.h"
#include "LinearSolvers/PETSc/petsc_solver.h"

#include "NeutronDiffusion/TransientSolver/transient_solver.h"

#include "Math/PETScUtils/petsc_utils.h"

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
  double magnitude = -0.01;
  double duration = 1.0;
  double interface = 40.0;

  std::string xsdir = "xs";
  std::string outdir = "outputs";

  for (int i = 0; i < argc; ++i)
  {
    std::string arg(argv[i]);
    std::cout << "Parsing argument " << i << " " << arg << std::endl;

    if (arg.find("magnitude") == 0)
      magnitude = std::stod(arg.substr(arg.find('=') + 1));
    else if (arg.find("duration") == 0)
      duration = std::stod(arg.substr(arg.find('=') + 1));
    else if (arg.find("interface") == 0)
      interface = std::stod(arg.substr(arg.find('=') + 1));
    else if (arg.find("output_directory") == 0)
      outdir = arg.substr(arg.find('=') + 1);
    else if (arg.find("xs_directory") == 0)
      xsdir = arg.substr(arg.find('=') + 1);
  }

  //============================================================
  // Mesh
  //============================================================

  std::vector<double> zones({0.0, interface, 200.0, 240.0});
  std::vector<size_t> n_cells({20, 80, 20});
  std::vector<unsigned int> materia_ids({0, 1, 2});
  auto mesh = create_1d_orthomesh(zones, n_cells, materia_ids);

  //============================================================
  // Materials
  //============================================================

  auto ramp_function =
      [magnitude, duration](const unsigned int group_num,
                            const std::vector<double>& args,
                            const double reference)
      {
        const double t = args[0];
        if (group_num == 1)
        {
          if (t > 0.0 && t <= duration)
            return (1.0 + t / duration * magnitude) * reference;
          else
            return (1.0 + magnitude) * reference;
        }
        else
          return reference;
      };

  std::vector<std::shared_ptr<Material>> materials;
  materials.emplace_back(std::make_shared<Material>("Material 0"));
  materials.emplace_back(std::make_shared<Material>("Material 1"));
  materials.emplace_back(std::make_shared<Material>("Material 2"));

  std::vector<std::shared_ptr<CrossSections>> xs;
  for (unsigned int i = 0; i < materials.size(); ++i)
    xs.emplace_back(std::make_shared<CrossSections>());

  std::vector<std::string> xs_paths;
  xs_paths.emplace_back(xsdir+"/fuel0.xs");
  xs_paths.emplace_back(xsdir+"/fuel1.xs");
  xs_paths.emplace_back(xsdir+"/fuel0.xs");

  for (unsigned int i = 0; i < materials.size(); ++i)
  {
    xs[i]->read_xs_file(xs_paths[i]);
    materials[i]->properties.emplace_back(xs[i]);
  }

  xs[0]->sigma_a_function = ramp_function;

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

  TransientSolver solver;

  solver.mesh = mesh;
  for (auto& material : materials)
    solver.materials.emplace_back(material);
  solver.linear_solver = linear_solver;

  solver.verbosity = 1;
  solver.use_precursors = true;

  solver.outer_tolerance = 1.0e-10;
  solver.max_outer_iterations = 1000;

  solver.algorithm = Algorithm::DIRECT;

  //============================================================
  // Define boundary conditions
  //============================================================

  solver.boundary_info.emplace_back(BoundaryType::ZERO_FLUX, -1);
  solver.boundary_info.emplace_back(BoundaryType::ZERO_FLUX, -1);

  //============================================================
  // Define transient parameters
  //============================================================

  solver.t_end = 2.0;
  solver.dt = 0.04;
  solver.time_stepping_method = TimeSteppingMethod::CRANK_NICHOLSON;

  solver.normalization_method = NormalizationMethod::TOTAL_POWER;
  solver.normalize_fission_xs = true;

  solver.write_outputs = true;
  solver.output_directory = outdir;

  solver.adaptive_time_stepping = true;
  solver.coarsen_threshold = 0.01;
  solver.refine_threshold = 0.05;

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

  std::cout << "\nSimulation Time: "
            << timer.get_time() << " ms\n";
}