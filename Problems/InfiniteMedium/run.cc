#include "Math/Quadratures/gauss_legendre.h"

#include "material.h"
#include "CrossSections/cross_sections.h"

#include "InfiniteMedium/KEigenvalueSolver/keigenvalue_solver.h"
#include "InfiniteMedium/TransientSolver/transient_solver.h"

#include <vector>
#include <cassert>

#include <iostream>
#include <fstream>
#include <filesystem>


using namespace PDEs;
using namespace Math;
using namespace Physics;
using namespace InfiniteMedium;


int main(int argc, char** argv)
{
  try{

    //##################################################
    // Parse the command line inputs
    //##################################################

    unsigned int n_grps = 618;
    std::string zaid = "1001";
    std::string ic = "maxwell_room";

    for (int i = 1; i < argc; ++i)
    {
      std::string arg(argv[i]);
      std::cout << "Parsing argument " << i << " " << arg << std::endl;

      if (arg.find("ng") == 0)
        n_grps = std::stoi(arg.substr(arg.find('=') + 1));
      else if (arg.find("zaid") == 0)
        zaid = arg.substr(arg.find('=') + 1);
      else if (arg.find("ic") == 0)
        ic = arg.substr(arg.find('=') + 1);
    }

    //##################################################
    // Format file paths
    //##################################################

    std::string filepath = __FILE__;
    std::string dirpath = filepath.substr(0, filepath.rfind('/') + 1);
    std::string xspath = dirpath + "xs/";
    std::string icpath = dirpath + "ics/";
    std::string outdir = dirpath + "outputs/";
    std::string grp = std::to_string(n_grps);


    // Path to cross sections
    xspath += zaid + "_" + grp + "g.ndi";
    assert(std::filesystem::is_regular_file(xspath));

    // Path to initial condition file
    icpath += grp + "g/" + ic + ".txt";
    assert(std::filesystem::is_regular_file(icpath));

    // Path to output directory
    outdir += grp + "g/" + ic;
    if (!std::filesystem::is_directory(outdir))
      std::filesystem::create_directories(outdir);

    //##################################################
    // Setup the problem
    //##################################################

    auto quadrature = std::make_shared<GaussLegendreQuadrature>(8);

    auto xs = std::make_shared<CrossSections>();
    xs->read_ndi_file(xspath);
    xs->make_pure_scatterer();

    std::vector<double> src_vals(xs->n_groups, 0.0);
    auto src = std::make_shared<IsotropicMultiGroupSource>(src_vals);

    TransientSolver solver;
    solver.xs = xs;
    solver.quadrature = quadrature;
    solver.src = src;

    solver.inner_tolerance = 1.0e-5;
    solver.max_inner_iterations = (unsigned int)1000;

    solver.verbosity = 2;

    // read in initial condition
    std::ifstream file(icpath);
    assert(file.is_open());

    std::string line;
    unsigned int group; double val;
    while (std::getline(file, line))
      if (line.find('#') != 0)
      {
        std::istringstream line_stream(line);
        line_stream >> group >> val;
        assert(group < xs->n_groups);
        solver.initial_conditions[group] = val;
      }
    file.close();

    solver.dt = 0.01;
    solver.t_end = 10 * solver.dt;

    solver.write_outputs = true;
    solver.output_directory = outdir;

    solver.initialize();
    solver.execute();
//
//    solver.write_angular_flux("Problems/InfiniteMedium/outputs");
//    solver.write_flux_moments("Problems/InfiniteMedium/outputs");
//    xs->write_group_structure("Problems/InfiniteMedium/outputs");
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