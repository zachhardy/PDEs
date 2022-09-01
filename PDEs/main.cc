#include "Math/Quadratures/gauss_legendre.h"

#include "material.h"
#include "CrossSections/cross_sections.h"

#include "InfiniteMedium/KEigenvalueSolver/keigenvalue_solver.h"
#include "InfiniteMedium/TransientSolver/transient_solver.h"

#include <iostream>
#include <vector>


using namespace PDEs;
using namespace Math;
using namespace Physics;
using namespace InfiniteMedium;


/**
 * The main execution function.
 * \param argc int      Number of supplied arguments.
 * \param argv char**   Array of strings for each argument.
 */
int main(int argc, char** argv)
{
  try{
    auto quadrature = std::make_shared<GaussLegendreQuadrature>(2);

    auto xs = std::make_shared<CrossSections>();
    xs->read_xs_file("xs/test_2g.xs");

    std::vector<double> src_vals(xs->n_groups, 0.0);
    auto src = std::make_shared<IsotropicMultiGroupSource>(src_vals);

    std::map<unsigned int, double> ics;
    ics[0] = 1.0;

    TransientSolver solver;
    solver.xs = xs;
    solver.quadrature = quadrature;
    solver.src = src;

    solver.initial_conditions = ics;

    solver.t_end = 100.0;
    solver.dt = 1.0;

    solver.inner_tolerance = 1.0e-6;
    solver.max_inner_iterations = (unsigned int)100;

//    solver.outer_tolerance = 1.0e-6;
//    solver.max_outer_iterations = 100;

    solver.use_dsa = false;
    solver.verbosity = 2;

    solver.initialize();
    solver.execute();
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