#include "mesh.h"
#include "ortho_grids.h"

#include "material.h"
#include "CrossSections/cross_sections.h"

#include "LinearSolvers/DirectSolvers"
#include "LinearSolvers/IterativeSolvers"
#include "LinearSolvers/PETSc/petsc_solver.h"

#include "NeutronDiffusion/Groupset/groupset.h"
#include "NeutronDiffusion/SteadyStateSolver/steadystate_solver.h"
#include "NeutronDiffusion/KEigenvalueSolver/keigenvalue_solver.h"

#include "macros.h"

#include <iostream>
#include <vector>

#include <petsc.h>

/**The main execution function.
 * \param argc int      Number of supplied arguments.
 * \param argv char**   Array of strings for each argument.
 */
int main(int argc, char** argv)
{
  try{
    using namespace Grid;

    using namespace Math;
    using namespace Math::LinearSolver;

    using namespace Physics;

    using namespace NeutronDiffusion;

    //================================================== Mesh

    size_t n_cells = 500;
    double slab_width = 1.0;
    double cell_width = slab_width / n_cells;

    std::vector<double> vertices(1, 0.0);
    for (size_t i = 0; i < n_cells; ++i)
      vertices.emplace_back(vertices.back() + cell_width);

    auto mesh = create_1d_mesh(vertices, CoordinateSystem::CARTESIAN);

    //================================================== Materials

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

    //================================================== Linear Solver

    Options opts;
    opts.verbosity = 1;
    opts.tolerance = 1.0e-6;

//    auto linear_solver = std::make_shared<PETScSolver>(KSPCG, PCNONE, opts);
    auto linear_solver = std::make_shared<CG>(opts);

    //================================================== Create the solver
    SteadyStateSolver solver;
    solver.mesh = mesh;
    solver.materials.emplace_back(material);
    solver.linear_solver = linear_solver;

    solver.verbose = true;
    solver.use_precursors = false;

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

    solver.solution_technique = SolutionTechnique::GROUPSET_WISE;


    //================================================== Run the problem

    PetscInitialize(&argc,&argv,(char*)0,NULL);

    solver.initialize();
    solver.execute();

    PetscFinalize();
  }
  catch (std::exception &exc) {
    std::cerr << std::endl
              << "----------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing:\n\n"
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
    std::cerr << "Unknown exception!\n\n"
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  return 0;
}