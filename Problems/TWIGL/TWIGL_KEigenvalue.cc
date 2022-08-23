#include "ortho_grids.h"

#include "material.h"
#include "CrossSections/cross_sections.h"

#include "timer.h"

#include "LinearSolvers/iterative_solvers.h"
#include "LinearSolvers/direct_solvers.h"
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
  double magnitude = 0.97667 - 1.0;
  double duration = 0.2;
  double sigs_01 = 0.01;
  std::string outdir = "Problems/TWIGL/outputs";

  for (int i = 0; i < argc; ++i)
  {
    std::string arg(argv[i]);
    std::cout << "Parsing argument " << i << " " << arg << std::endl;

    if (arg.find("magnitude") == 0)
      magnitude = std::stod(arg.substr(arg.find("=") + 1));
    else if (arg.find("duration") == 0)
      duration = std::stod(arg.substr(arg.find("=") + 1));
    else if (arg.find("scatter") == 0)
      sigs_01 = std::stod(arg.substr(arg.find("=") + 1));
    else if (arg.find("output_directory") == 0)
      outdir = arg.substr(arg.find("=") + 1);
  }

  //============================================================
  // Mesh
  //============================================================

  size_t n_x = 41, n_y = 41;
  double X = 80.0, Y = 80.0;
  double dx = X / (n_x - 1), dy = Y / (n_y - 1);

  std::vector<double> x_verts(1, 0.0);
  for (size_t i = 0; i < n_x - 1; ++i)
    x_verts.push_back(x_verts.back() + dx);

  std::vector<double> y_verts(1, 0.0);
  for (size_t i = 0; i < n_y - 1; ++i)
    y_verts.push_back(y_verts.back() + dy);

  auto mesh = create_2d_orthomesh(x_verts, y_verts);

  for (auto& cell : mesh->cells)
  {
    auto& c = cell.centroid;
    if ((c.x() > 24.0 and c.x() < 56.0) and (c.y() > 24.0 and c.y() < 56.0))
      cell.material_id = 0;
    else if (c.x() < 24.0 and (c.y() > 24.0 and c.y() < 56.0) or
             ((c.x() > 24.0 and c.x() < 56.0) and c.y() < 24.0))
      cell.material_id = 1;
    else
      cell.material_id = 2;
  }

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
        if (t >= 0.0 and t <= duration)
          return (1.0 + t / duration * magnitude) * reference;
        else
          return (1.0 + magnitude) * reference;
      }
      else
        return reference;
    };

  std::vector<std::shared_ptr<Material>> materials;
  materials.emplace_back(std::make_shared<Material>("Fuel 0 Perturbed"));
  materials.emplace_back(std::make_shared<Material>("Fuel 0"));
  materials.emplace_back(std::make_shared<Material>("Fuel 1"));

  std::vector<std::shared_ptr<CrossSections>> xs;
  for (unsigned int i = 0; i < materials.size(); ++i)
    xs.emplace_back(std::make_shared<CrossSections>());

  std::vector<std::string> xs_paths;
  xs_paths.emplace_back("Problems/TWIGL/xs/fuel0.xs");
  xs_paths.emplace_back("Problems/TWIGL/xs/fuel0.xs");
  xs_paths.emplace_back("Problems/TWIGL/xs/fuel1.xs");

  for (unsigned int i = 0; i < materials.size(); ++i)
  {
    xs[i]->read_xs_file(xs_paths[i]);
    if (i == 2)
      xs[i]->transfer_matrices[0][1][0] = sigs_01;
    materials[i]->properties.emplace_back(xs[i]);
  }

  xs[0]->sigma_a_function = ramp_function;

  const unsigned int n_groups = xs.front()->n_groups;

  //============================================================
  // Linear Solver
  //============================================================

  Options opts;
  opts.verbosity = 0;
  opts.tolerance = 1.0e-14;
  opts.max_iterations = 10000;

  std::shared_ptr<LinearSolverBase<SparseMatrix>> linear_solver;
  linear_solver = std::make_shared<PETScSolver>(KSPCG, PCLU, opts);

  //============================================================
  // Create the diffusion solver
  //============================================================

  KEigenvalueSolver solver;

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
  // Initialize groups and groupsets
  //============================================================

  for (unsigned int g = 0; g < n_groups; ++g)
    solver.groups.emplace_back(g);

  //============================================================
  // Define boundary conditions
  //============================================================

  solver.boundary_info.emplace_back(BoundaryType::REFLECTIVE, -1);
  solver.boundary_info.emplace_back(BoundaryType::VACUUM, -1);
  solver.boundary_info.emplace_back(BoundaryType::VACUUM, -1);
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
