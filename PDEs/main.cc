#include "grid_structs.h"

#include "material.h"
#include "CrossSections/cross_sections.h"

#include "math.h"
#include "linear_solver.h"

#include "NeutronDiffusion/Groupset/groupset.h"
#include "NeutronDiffusion/SteadyStateSolver/FV/steadystate_solver_fv.h"

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


    //================================================== Create the mesh
    size_t n_cells = 100;
    double slab_width = 1.0;
    double cell_width = slab_width / n_cells;

    std::vector<double> vertices(1, 0.0);
    for (size_t i = 0; i < n_cells; ++i)
      vertices.emplace_back(vertices.back() + cell_width);

    auto mesh = create_1d_mesh(vertices, CoordinateSystem::CARTESIAN);

    //================================================== Create the materials
    auto material = std::make_shared<Material>();

    // Create the cross sections
    auto xs = std::make_shared<CrossSections>();
    xs->read_xs_file("xs_data/test_1g.xs");
    material->properties.emplace_back(xs);

    size_t n_groups = xs->n_groups;

    // Create the multigroup source
    std::vector<double> mg_source(n_groups, 1.0);
    auto src = std::make_shared<IsotropicMultiGroupSource>(mg_source);
    material->properties.emplace_back(src);

    //================================================== Create the solver
    SteadyStateSolver_FV solver;
    solver.mesh = mesh;
    solver.materials.emplace_back(material);

    // Initialize groups
    for (size_t g = 0; g < n_groups; ++g)
      solver.groups.emplace_back(g);

    // Initialize groupsets
    Groupset groupset(0);
    for (size_t g = 0; g < n_groups; ++g)
      groupset.groups.emplace_back(solver.groups[g]);
    solver.groupsets.emplace_back(groupset);

    // Define boundary conditions
    solver.boundary_info.emplace_back(BoundaryType::ZERO_FLUX, -1);
    solver.boundary_info.emplace_back(BoundaryType::ZERO_FLUX, -1);

    solver.options.solution_technique = SolutionTechnique::GROUPSET_WISE;


    //================================================== Run the problem
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