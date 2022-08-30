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


int main()
{
  try{
    auto quadrature = std::make_shared<GaussLegendreQuadrature>(8);

    auto xs = std::make_shared<CrossSections>();
    xs->read_ndi_file("xs_data/1001.ndi");

    std::vector<double> src_vals(xs->n_groups, 0.0);
    src_vals[xs->n_groups - 2] = 1.0;
    auto src = std::make_shared<IsotropicMultiGroupSource>(src_vals);


    SteadyStateSolver solver;
    solver.xs = xs;
    solver.quadrature = quadrature;
    solver.src = src;

    solver.inner_tolerance = 1.0e-6;
    solver.max_inner_iterations = (unsigned int)10000;

    solver.use_dsa = false;
    solver.verbosity = 2;

//    std::map<unsigned int, double> ics;
//    for (unsigned int g = 0; g < xs->n_groups; ++g)
//      solver.inital_conditions[g] = 1.0;
//
//    solver.t_end = 100.0;
//    solver.dt = 1.0;

    solver.initialize();
    solver.execute();

    solver.write_angular_flux("psi");
    solver.write_flux_moments("phi");
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