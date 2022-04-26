#include "grid_structs.h"
#include "Discretization/FiniteVolume/finite_volume.h"

#include "material.h"
#include "CrossSections/cross_sections.h"

#include "math.h"
#include "linear_solver.h"

#include "NeutronDiffusion/SteadyStateSolver/steadystate_solver.h"

#include <iostream>
#include <vector>

/**The main execution function.
 * \param argc int      Number of supplied arguments.
 * \param argv char**   Array of strings for each argument.
 */
int main(int argc, char** argv)
{
  try{
    using namespace grid;
    using namespace math;
    using namespace physics;
    using namespace neutron_diffusion;

    // Create the solver
    SteadyStateSolver solver;

    // Create the vertices
    size_t n_cells = 100;
    double slab_width = 1.0;
    double cell_width = slab_width / n_cells;

    std::vector<double> vertices(1, 0.0);
    for (size_t i = 0; i < n_cells; ++i)
      vertices.emplace_back(vertices.back() + cell_width);

    auto mesh = create_1d_mesh(vertices, CoordinateSystem::CARTESIAN);
    solver.mesh = mesh;

    auto discretization = std::make_shared<FiniteVolume>(mesh);
    solver.discretization = discretization;

    auto material = std::make_shared<Material>();

    auto xs = std::make_shared<CrossSections>();
    xs->read_xs_file("xs_data/test_2g.xs");
    material->properties.emplace_back(xs);

    std::vector<double> mg_source(xs->n_groups, 1.0);
    auto src = std::make_shared<IsotropicMultiGroupSource>(mg_source);
    material->properties.emplace_back(src);

    solver.materials.emplace_back(material);

    solver.boundary_info.emplace_back(BoundaryType::ZERO_FLUX, -1);
    solver.boundary_info.emplace_back(BoundaryType::ZERO_FLUX, -1);

    solver.solution_method = SolutionMethod::ITERATIVE;
    solver.linear_solver_type = LinearSolverType::LU;


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