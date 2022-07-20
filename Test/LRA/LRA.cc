#include "mesh.h"
#include "ortho_grids.h"

#include "material.h"
#include "CrossSections/cross_sections.h"

#include "timer.h"

#include "LinearSolvers/DirectSolvers"
#include "LinearSolvers/IterativeSolvers"
#include "LinearSolvers/PETSc/petsc_solver.h"

#include "NeutronDiffusion/groupset.h"
#include "NeutronDiffusion/TransientSolver/transient_solver.h"

#include <iostream>
#include <vector>

#include <petsc.h>


int main(int argc, char** argv)
{
  //============================================================
  // Mesh
  //============================================================
  using namespace Grid;

  size_t n_x = 23, n_y = 23;
  double X = 165.0, Y = 165.0;
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
    if (c.x < 15.0)
    {
      if (c.y < 15.0) cell.material_id = 1;
      else if (c.y > 15.0 and c.y < 75.0) cell.material_id = 0;
      else if (c.y > 75.0 and c.y < 105.0) cell.material_id = 1;
      else if (c.y > 105.0 and c.y < 135.0) cell.material_id = 2;
      else cell.material_id = 5;
    }
    else if (c.x > 15.0 and c.x < 75.0)
    {
      if (c.y < 105.0) cell.material_id = 0;
      else if (c.y > 105.0 and c.y < 135.0) cell.material_id = 2;
      else cell.material_id = 5;
    }
    else if (c.x > 75.0 and c.x < 105.0)
    {
      if (c.y < 15.0) cell.material_id = 1;
      else if (c.y > 15.0 and c.y < 75.0) cell.material_id = 0;
      else if (c.y > 75.0 and c.y < 105.0) cell.material_id = 1;
      else if (c.y > 105.0 and c.y < 135.0) cell.material_id = 2;
      else cell.material_id = 5;
    }
    else if (c.x > 105.0 and c.x < 120.0)
    {
      if (c.y < 75.0) cell.material_id = 2;
      else if (c.y > 75.0 and c.y < 105.0) cell.material_id = 4;
      else if (c.y > 105.0 and c.y < 120.0) cell.material_id = 3;
      else cell.material_id = 5;
    }
    else if (c.x > 120.0 and c.x < 135.0)
    {
      if (c.y < 75.0) cell.material_id = 2;
      else if (c.y > 75.0 and c.y < 105.0) cell.material_id = 4;
      else cell.material_id = 5;
    }
    else
      cell.material_id = 5;
  }

  //============================================================
  // Materials
  //============================================================
  using namespace Physics;

  const double delta = 0.8787631 - 1.0;
  const double t_ramp = 2.0;
  const double gamma = 3.034e-3;

  auto rod_ejection_with_feedback =
      [delta, t_ramp, gamma](
          const unsigned int group_num,
          const double current_time,
          const double temperature,
          const double reference_temperature,
          const double reference_value)
      {
        if (group_num == 0)
        {
          const double sqrt_T = std::sqrt(temperature);
          const double sqrt_T0 = std::sqrt(reference_temperature);
          return reference_value*(1.0 + gamma*(sqrt_T - sqrt_T0));
        }
        else if (group_num == 1)
        {
          if (current_time <= t_ramp)
            return (1.0 + current_time/t_ramp * delta) * reference_value;
          else
            return (1.0 + delta) * reference_value;
        }
        else
          return reference_value;
      };

  auto feedback =
      [gamma](
          const unsigned int group_num,
          const double current_time,
          const double temperature,
          const double reference_temperature,
          const double reference_value)
      {
        if (group_num == 0)
        {
          const double sqrt_T = std::sqrt(temperature);
          const double sqrt_T0 = std::sqrt(reference_temperature);
          return reference_value*(1.0 + gamma*(sqrt_T - sqrt_T0));
        }
        else
          return reference_value;
      };


  std::vector<std::shared_ptr<Material>> materials;
  materials.emplace_back(std::make_shared<Material>("Fuel 1 w/ Rod"));
  materials.emplace_back(std::make_shared<Material>("Fuel 1 w/o Rod"));
  materials.emplace_back(std::make_shared<Material>("Fuel 2 w/ Rod"));
  materials.emplace_back(std::make_shared<Material>("Fuel 2 w/o Rod"));
  materials.emplace_back(std::make_shared<Material>("Rod Ejection Region"));
  materials.emplace_back(std::make_shared<Material>("Reflector"));

  std::vector<std::shared_ptr<CrossSections>> xs;
  for (size_t i = 0; i < materials.size(); ++i)
    xs.emplace_back(std::make_shared<CrossSections>());

  std::vector<std::string> xs_paths;
  xs_paths.emplace_back("Test/LRA/xs/fuel1_w_rod.xs");
  xs_paths.emplace_back("Test/LRA/xs/fuel1_wo_rod.xs");
  xs_paths.emplace_back("Test/LRA/xs/fuel2_w_rod.xs");
  xs_paths.emplace_back("Test/LRA/xs/fuel2_wo_rod.xs");
  xs_paths.emplace_back("Test/LRA/xs/fuel2_w_rod.xs");
  xs_paths.emplace_back("Test/LRA/xs/reflector.xs");

  for (size_t i = 0; i < materials.size(); ++i)
  {
    xs[i]->read_xs_file(xs_paths[i]);

    if (i < 4) xs[i]->sigma_a_function = feedback;
    else if (i == 4) xs[i]->sigma_a_function = rod_ejection_with_feedback;

    materials[i]->properties.emplace_back(xs[i]);
  }

  const size_t n_groups = xs.front()->n_groups;

  //============================================================
  // Linear Solver
  //============================================================
  using namespace LinearSolver;

  Options opts;
  opts.verbosity = 0;
  opts.tolerance = 1.0e-14;
  opts.max_iterations = 10000;

  std::shared_ptr<LinearSolverBase<SparseMatrix>> linear_solver;
  linear_solver = std::make_shared<PETScSolver>(KSPGMRES, PCLU, opts);

  //============================================================
  // Create the diffusion solver
  //============================================================
  using namespace NeutronDiffusion;
  TransientSolver solver;
  solver.mesh = mesh;

  for (auto& material : materials)
    solver.materials.emplace_back(material);

  solver.linear_solver = linear_solver;

  solver.verbosity = 0;
  solver.use_precursors = true;

  solver.tolerance = 1.0e-12;
  solver.max_iterations = 1000;

  solver.solution_technique = SolutionTechnique::FULL_SYSTEM;

  //============================================================
  // Define transient parameters
  //============================================================
  solver.t_end = 3.0;
  solver.dt = 0.01;
  solver.time_stepping_method = TimeSteppingMethod::CRANK_NICHOLSON;

  solver.normalize_fission_xs = true;
  solver.normalization_method = NormalizationMethod::AVERAGE_POWER;
  solver.power = 1.0e-6;

  solver.write_outputs = true;
  solver.output_directory = "Test/LRA/outputs";

  solver.adaptivity = true;
  solver.coarsen_threshold = 0.01;
  solver.refine_threshold = 0.1;

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
