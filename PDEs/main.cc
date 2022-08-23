#include "Math/Polynomials/polynomial.h"

#include <iostream>
#include <fstream>
#include <vector>


/**The main execution function.
 * \param argc int      Number of supplied arguments.
 * \param argv char**   Array of strings for each argument.
 */
int main(int argc, char** argv)
{
  try{
    using namespace PDEs::Math::Polynomials;

    size_t n_pts = 101;
    unsigned int max_order = 3;

    std::vector<double> x(n_pts);
    x[0] = -1.0;
    for (size_t i = 1; i < n_pts - 1; ++i)
      x[i] = x[i - 1] + 2.0/(n_pts - 1);
    x[n_pts - 1] = 1.0;

    std::vector<std::vector<double>> y(max_order + 1);
    for (unsigned int p = 0; p <= max_order; ++p)
      for (size_t i = 0; i < n_pts; ++i)
        y[p].push_back(legendre(p, x[i]));

    std::ofstream ofile;
    ofile.open("plot.py");

    ofile << "import numpy as np\n"
          << "import matplotlib.pyplot as plt\n\n";

    ofile << "x = np.zeros(" << n_pts << ")\n";
    for (size_t i = 0; i < n_pts; ++i)
      ofile << "x[" << i << "] = " << x[i] << "\n";

    ofile << "y = np.zeros((" << max_order + 1 << ", " << n_pts << "))\n";
    for (unsigned int p = 0; p <= max_order; ++p)
      for (size_t i = 0; i < n_pts; ++i)
        ofile << "y[" << p << "][" << i << "] = " << y[p][i] << "\n";
    ofile << "\n";

    ofile << "plt.figure()\n"
          << "for p in range(" << max_order + 1 << "):\n"
          << "\tplt.plot(x, y[p], label=f'L={p}')\n"
          << "plt.legend()\n"
          << "plt.show()\n";

    ofile.close();
    system("python plot.py");

    auto leg_roots = legendre_roots(32);
    auto cheb_roots = chebyshev_roots(32);

    for (unsigned int i = 0; i < leg_roots.size(); ++i)
      std::cout << leg_roots[i] << "  " << cheb_roots[i] << std::endl;
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