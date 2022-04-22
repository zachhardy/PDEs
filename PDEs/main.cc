#include "grid_structs.h"
#include "FiniteVolume/finite_volume.h"

#include "material.h"
#include "CrossSections/cross_sections.h"

#include "math.h"
#include "linear_solver.h"

#include "Diffusion/SteadyStateSolver/steadystate_solver.h"

#include <iostream>
#include <vector>

/**The main execution function.
 * \param argc int      Number of supplied arguments.
 * \param argv char**   Array of strings for each argument.
 */
int main(int argc, char** argv)
{
  try{
    using namespace diffusion;

    // Create solver
    auto solver = SteadyStateSolver();

    // Create the mesh
    auto coord_sys = grid::CoordinateSystem::CARTESIAN;
    std::vector<double> zone_edges = {0.0, 1.0};
    std::vector<size_t> zone_subdivisions = {4};
    std::vector<int> material_ids = {0};

    // Create mesh
    solver.mesh = create_1d_mesh(zone_edges, zone_subdivisions,
                                 material_ids, coord_sys);

    // Create discretization
    typedef discretization::FiniteVolume FV;
    solver.discretization = std::make_shared<FV>(solver.mesh);

    // Create materials
    solver.materials.push_back(std::make_shared<material::Material>());

    // Cross sections
    auto xs = std::make_shared<material::CrossSections>();
    xs->read_xs_file("xs_data/test.xs");
    solver.materials.back()->properties.push_back(xs);

    // Multigroup source
    typedef material::IsotropicMultiGroupSource IsotropicMultiGroupSource;
    std::vector<double> mg_source = {1.0};
    auto src = std::make_shared<IsotropicMultiGroupSource>(mg_source);
    solver.materials.back()->properties.push_back(src);

    // Boundaries
    std::vector<std::vector<double>> bndry_vals;
    solver.boundary_info.emplace_back(BoundaryType::ZERO_FLUX, bndry_vals);
    solver.boundary_info.emplace_back(BoundaryType::ZERO_FLUX, bndry_vals);

    // Linear solver option
    solver.linear_solver_type = linear_solver::LinearSolverType::LU;

    // Run the problem
    solver.initialize();
    solver.execute();
  }
  catch (std::exception &exc) {
    std::cerr << std::endl
              << "----------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------"
              << std::endl;
    return 1;
  }
  catch (...) {
    std::cerr << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  return 0;
}