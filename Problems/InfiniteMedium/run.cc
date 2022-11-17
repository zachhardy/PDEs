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


using namespace std;
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
    string zaid = "1001";
    string ic = "maxwell_room";

    for (int i = 1; i < argc; ++i)
    {
      string arg(argv[i]);
      cout << "Parsing argument " << i << " " << arg << endl;

      if (arg.find("ng") == 0)
        n_grps = stoi(arg.substr(arg.find('=') + 1));
      else if (arg.find("zaid") == 0)
        zaid = arg.substr(arg.find('=') + 1);
      else if (arg.find("ic") == 0)
        ic = arg.substr(arg.find('=') + 1);
    }

    //##################################################
    // Format file paths
    //##################################################

    string filepath = __FILE__;
    string dirpath = filepath.substr(0, filepath.rfind('/') + 1);
    string xspath = dirpath + "xs/";
    string icpath = dirpath + "ics/";
    string outdir = dirpath + "outputs/";
    string grp = to_string(n_grps);


    // Path to cross sections
    xspath += zaid + "_" + grp + "g.ndi";
    assert(filesystem::is_regular_file(xspath));

    // Path to initial condition file
    icpath += grp + "g/" + ic + ".txt";
    assert(filesystem::is_regular_file(icpath));

    // Path to output directory
    outdir += grp + "g/" + ic;
    if (!filesystem::is_directory(outdir))
      filesystem::create_directories(outdir);

    //##################################################
    // Setup the problem
    //##################################################

    auto quadrature = make_shared<GaussLegendreQuadrature>(8);

    auto xs = make_shared<CrossSections>();
    xs->read_ndi_file(xspath);
    xs->make_pure_scatterer();

    vector<double> src_vals(xs->n_groups, 0.0);
    auto src = make_shared<IsotropicMultiGroupSource>(src_vals);

    TransientSolver solver;
    solver.xs = xs;
    solver.quadrature = quadrature;
    solver.src = src;

    solver.inner_tolerance = 1.0e-5;
    solver.max_inner_iterations = (unsigned int)1000;

    solver.verbosity = 2;

    // read in initial condition
    ifstream file(icpath);
    assert(file.is_open());

    string line;
    unsigned int group; double val;
    while (getline(file, line))
      if (line.find('#') != 0)
      {
        istringstream line_stream(line);
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
  catch (exception &exc) {
    cerr << endl
              << "----------------------------------------"
              << endl;
    cerr << "Exception on processing:\n\n"
              << exc.what() << endl
              << "Aborting!" << endl
              << "----------------------------------------"
              << endl;
    return 1;
  }
  catch (...) {
    cerr << endl
              << "----------------------------------------------------"
              << endl;
    cerr << "Unknown exception!\n\n"
              << "Aborting!" << endl
              << "----------------------------------------------------"
              << endl;
    return 1;
  }
  return 0;
}