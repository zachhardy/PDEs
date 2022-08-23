#include "Math/Polynomials/polynomial.h"
#include "Math/Quadratures/gauss_legendre.h"

#include "material.h"
#include "CrossSections/cross_sections.h"

#include "NeutronTransport/InfiniteMedium/infinite_medium.h"

#include <iostream>
#include <fstream>
#include <vector>


void
generate_grid(const size_t n_pts, std::vector<double>& x)
{
  x.resize(n_pts, 0.0);
  x[0] = -1.0;
  for (size_t i = 1; i < n_pts - 1; ++i)
    x[i] = x[i - 1] + 2.0/(n_pts - 1);
  x[n_pts - 1] = 1.0;
}


void
eval_legendre(const unsigned int max_order,
              const std::vector<double>& x,
              std::vector<std::vector<double>>& y)
{
  using namespace PDEs::Math::Polynomials;
  y.resize(max_order + 1, std::vector<double>(x.size()));
  for (unsigned int p = 0; p <= max_order; ++p)
    for (size_t i = 0; i < x.size(); ++i)
      y[p][i] = legendre(p, x[i]);
}


void
eval_chebyshev(const unsigned int max_order,
               const std::vector<double>& x,
               std::vector<std::vector<double>>& y)
{
  using namespace PDEs::Math::Polynomials;
  y.resize(max_order + 1, std::vector<double>(x.size()));
  for (unsigned int p = 0; p <= max_order; ++p)
    for (size_t i = 0; i < x.size(); ++i)
      y[p][i] = chebyshev(p, x[i]);
}


void
plot(const std::vector<double>& x,
     const std::vector<std::vector<double>>& y,
     const std::string title = "")
{
  std::ofstream ofile;
  ofile.open("plot.py");

  ofile << "import numpy as np\n"
        << "import matplotlib.pyplot as plt\n\n";

  ofile << "x = np.zeros(" << x.size() << ")\n";
  for (size_t i = 0; i < x.size(); ++i)
    ofile << "x[" << i << "] = " << x[i] << "\n";

  ofile << "y = np.zeros((" << y.size() << ", " << x.size() << "))\n";
  for (unsigned int p = 0; p < y.size(); ++p)
    for (size_t i = 0; i < x.size(); ++i)
      ofile << "y[" << p << "][" << i << "] = " << y[p][i] << "\n";
  ofile << "\n";

  ofile << "plt.figure()\n"
        << "plt.title(\"" << title << "\")\n"
        << "for p in range(" << y.size() << "):\n"
        << "\tplt.plot(x, y[p], label=f'L={p}')\n"
        << "plt.legend()\n"
        << "plt.show()\n";

  ofile.close();
  system("python plot.py");
}


void
plot_legendre(const size_t n_pts = 101, const unsigned int max_order = 3)
{
  using namespace PDEs::Math::Polynomials;

  std::vector<double> x;
  generate_grid(n_pts, x);

  std::vector<std::vector<double>> y;
  eval_legendre(max_order, x, y);

  plot(x, y, "Legendre Polynomials");
}


void
plot_chebyshev(const size_t n_pts = 100, const unsigned int max_order = 3)
{
  using namespace PDEs::Math::Polynomials;

  std::vector<double> x;
  generate_grid(n_pts, x);

  std::vector<std::vector<double>> y;
  eval_chebyshev(max_order, x, y);

  plot(x, y, "Chebyshev Polynomials");
}


using namespace PDEs;
using namespace Math;
using namespace NeutronTransport;


/**The main execution function.
 * \param argc int      Number of supplied arguments.
 * \param argv char**   Array of strings for each argument.
 */
int main(int argc, char** argv)
{
  try{

    auto quadrature = std::make_shared<GaussLegendreQuadrature>(8);

    auto xs = std::make_shared<CrossSections>();
    xs->read_xs_file("xs_data/graphite_pure.xs");

    std::vector<double> src_vals(xs->n_groups, 0.0);
    src_vals[xs->n_groups - 1] = 1.0;
    auto src = std::make_shared<IsotropicMultiGroupSource>(src_vals);

    InfiniteMedium solver;
    solver.xs = xs;
    solver.quadrature = quadrature;
    solver.src = src;

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