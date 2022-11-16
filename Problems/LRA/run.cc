#include "ortho_grids.h"

#include "material.h"
#include "CrossSections/cross_sections.h"

#include "timer.h"

#include "Math/LinearSolvers/Iterative/cg.h"
#include "Math/LinearSolvers/Direct/cholesky.h"
#include "LinearSolvers/PETSc/petsc_solver.h"

#include "NeutronDiffusion/TransientSolver/transient_solver.h"

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
  double magnitude = 0.8787631 - 1.0;
  double duration = 2.0;
  double feedback = 3.034e-3;

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
    else if (arg.find("feedback") == 0)
      feedback = std::stod(arg.substr(arg.find('=') + 1));
    else if (arg.find("output_directory") == 0)
      outdir = arg.substr(arg.find('=') + 1);
    else if (arg.find("xs_directory") == 0)
      xsdir = arg.substr(arg.find('=') + 1);
  }

  //============================================================
  // Mesh
  //============================================================

  size_t n_x = 23, n_y = 23;
  double X = 165.0, Y = 165.0;
  double dx = X / (n_x - 1), dy = Y / (n_y - 1);

  std::vector<double> x_verts(1, 0.0);
  for (size_t i = 0; i < n_x - 1; ++i)
    x_verts.push_back(x_verts.back() + dx);

  std::vector<double> y_verts(1, 0.0);
  for (size_t i = 0; i < n_y - 1; ++i)
    y_verts.push_back(y_verts.back() + dy);

  auto mesh = create_2d_orthomesh(x_verts, y_verts);

  for (auto& cell: mesh->cells)
  {
    auto& c = cell.centroid;
    if (c.x() < 15.0)
    {
      if (c.y() < 15.0) cell.material_id = 1;
      else if (c.y() > 15.0 and c.y() < 75.0) cell.material_id = 0;
      else if (c.y() > 75.0 and c.y() < 105.0) cell.material_id = 1;
      else if (c.y() > 105.0 and c.y() < 135.0) cell.material_id = 2;
      else cell.material_id = 5;
    } else if (c.x() > 15.0 and c.x() < 75.0)
    {
      if (c.y() < 105.0) cell.material_id = 0;
      else if (c.y() > 105.0 and c.y() < 135.0) cell.material_id = 2;
      else cell.material_id = 5;
    } else if (c.x() > 75.0 and c.x() < 105.0)
    {
      if (c.y() < 15.0) cell.material_id = 1;
      else if (c.y() > 15.0 and c.y() < 75.0) cell.material_id = 0;
      else if (c.y() > 75.0 and c.y() < 105.0) cell.material_id = 1;
      else if (c.y() > 105.0 and c.y() < 135.0) cell.material_id = 2;
      else cell.material_id = 5;
    } else if (c.x() > 105.0 and c.x() < 120.0)
    {
      if (c.y() < 75.0) cell.material_id = 2;
      else if (c.y() > 75.0 and c.y() < 105.0) cell.material_id = 4;
      else if (c.y() > 105.0 and c.y() < 120.0) cell.material_id = 3;
      else cell.material_id = 5;
    } else if (c.x() > 120.0 and c.x() < 135.0)
    {
      if (c.y() < 75.0) cell.material_id = 2;
      else if (c.y() > 75.0 and c.y() < 105.0) cell.material_id = 4;
      else cell.material_id = 5;
    } else
      cell.material_id = 5;
  }

  //============================================================
  // Materials
  //============================================================

  const double T0 = 300.0;

  auto rod_ejection_with_feedback =
      [magnitude, duration, feedback, T0](
          const unsigned int group_num,
          const std::vector<double>& args,
          const double reference)
      {
        const double t = args[0], T = args[1];

        if (group_num == 0)
          return (1.0 + feedback * (std::sqrt(T) - std::sqrt(T0))) * reference;
        else if (group_num == 1)
        {
          if (t <= duration)
            return (1.0 + t / duration * magnitude) * reference;
          else
            return (1.0 + magnitude) * reference;
        }
        else
          return reference;
      };

  auto feedback_function =
      [feedback, T0](
          const unsigned int group_num,
          const std::vector<double>& args,
          const double reference)
      {
        const double T = args[1];

        if (group_num == 0)
          return (1.0 + feedback * (std::sqrt(T) - std::sqrt(T0))) * reference;
        else
          return reference;
      };


  std::vector<std::shared_ptr<Material>> materials;
  materials.emplace_back(std::make_shared<Material>("Fuel 1 w/ Rod"));
  materials.emplace_back(std::make_shared<Material>("Fuel 1 w/o Rod"));
  materials.emplace_back(std::make_shared<Material>("Fuel 2 w/ Rod"));
  materials.emplace_back(std::make_shared<Material>("Fuel 2 w/o Rod"));
  materials.emplace_back(std::make_shared<Material>("Rod Ejection Region"));
  materials.emplace_back(std::make_shared<Material>("Reflector"));

  std::vector<std::shared_ptr<CrossSections>> xs;
  for (unsigned int i = 0; i < materials.size(); ++i)
    xs.emplace_back(std::make_shared<CrossSections>());

  std::vector<std::string> xs_paths;
  xs_paths.emplace_back(xsdir+"/fuel1_w_rod.xs");
  xs_paths.emplace_back(xsdir+"/fuel1_wo_rod.xs");
  xs_paths.emplace_back(xsdir+"/fuel2_w_rod.xs");
  xs_paths.emplace_back(xsdir+"/fuel2_wo_rod.xs");
  xs_paths.emplace_back(xsdir+"/fuel2_w_rod.xs");
  xs_paths.emplace_back(xsdir+"/reflector.xs");

  for (unsigned int i = 0; i < materials.size(); ++i)
  {
    xs[i]->read_xs_file(xs_paths[i]);
    materials[i]->properties.emplace_back(xs[i]);
    if (i < 4) xs[i]->sigma_a_function = feedback_function;
    else if (i == 4) xs[i]->sigma_a_function = rod_ejection_with_feedback;
  }

  //============================================================
  // Linear Solver
  //============================================================

  Options opts;
  opts.verbosity = 0;
  opts.tolerance = 1.0e-14;
  opts.max_iterations = 10000;

  std::shared_ptr<LinearSolverBase<SparseMatrix>> linear_solver;
  linear_solver = std::make_shared<PETScSolver>(KSPGMRES, PCLU, opts);

  //============================================================
  // Create the diffusion solver
  //============================================================

  TransientSolver solver;
  solver.mesh = mesh;

  for (auto& material: materials)
    solver.materials.emplace_back(material);

  solver.linear_solver = linear_solver;

  solver.verbosity = 0;
  solver.use_precursors = true;

  solver.outer_tolerance = 1.0e-12;
  solver.max_outer_iterations = 1000;

  solver.algorithm = Algorithm::DIRECT;

  //============================================================
  // Define transient parameters
  //============================================================

  solver.t_end = 3.0;
  solver.dt = 0.01;
  solver.time_stepping_method = TimeSteppingMethod::CRANK_NICHOLSON;

  solver.normalize_fission_xs = true;
  solver.normalization_method = NormalizationMethod::AVERAGE_POWER;
  solver.initial_power = 1.0e-6;

  solver.write_outputs = true;
  solver.output_directory = outdir;

  solver.adaptive_time_stepping = true;
  solver.coarsen_threshold = 0.01;
  solver.refine_threshold = 0.1;

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
