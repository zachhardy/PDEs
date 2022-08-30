#include "cross_sections.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <numeric>
#include <cassert>


using namespace PDEs;
using namespace Physics;


void
CrossSections::
read_xs_file(const std::string file_name,
             const double rho,
             const bool verbose)
{
  std::cout << "Reading cross-section file \"" << file_name << "\".\n";

  // Clear old cross sections
  reset();

  density = rho;

  // Open the file
  std::ifstream file;
  file.open(file_name);
  assert(file.is_open());

  //############################################################
  /** Lambda for reading vector data. */
  auto read_1d_data =
      [](const std::string keyword,
         std::vector<double>& destination,
         const unsigned int n,
         std::ifstream& file,
         std::istringstream& line_stream,
         size_t& line_number)
      {
        std::string line;
        std::getline(file, line);
        line_stream = std::istringstream(line);
        ++line_number;

        unsigned int count = 0;
        int i; double value;
        while (line != keyword + "_END")
        {
          line_stream >> i >> value;
          destination.at(i) = value;
          assert(count++ <= n);

          std::getline(file, line);
          line_stream = std::istringstream(line);
          ++line_number;
        }
      };

  //############################################################
  /** Lambda for reading transfer matrix data. */
  auto read_transfer_matrix =
      [](const std::string keyword,
         std::vector<TransferMatrix>& destination,
         std::ifstream& file,
         std::istringstream& line_stream,
         size_t& line_number)
      {
        std::string word, line;
        std::getline(file, line);
        line_stream = std::istringstream(line);
        ++line_number;

        int ell, g, gp;
        double value;
        while (line != keyword + "_END")
        {
          line_stream >> word;
          if (word == "M_GPRIME_G_VAL")
          {
            line_stream >> ell >> gp >> g >> value;
            destination.at(ell).at(g).at(gp) = value;
          }

          std::getline(file, line);
          line_stream = std::istringstream(line);
          ++line_number;
        }
      };

  //############################################################
  /** Lambda for reading delayed emission spectra data. */
  auto read_chi_delayed =
      [](const std::string keyword,
         EmissionSpectra& destination,
         std::ifstream& file,
         std::istringstream& line_stream,
         size_t& line_number)
      {
        std::string word, line;
        std::getline(file, line);
        line_stream = std::istringstream(line);
        ++line_number;

        int g, j;
        double value;
        while (line != keyword + "_END")
        {
          line_stream >> word;
          if (word == "G_PRECURSORJ_VAL")
          {
            line_stream >> g >> j >> value;
            destination.at(g).at(j) = value;
          }

          std::getline(file, line);
          line_stream = std::istringstream(line);
          ++line_number;
        }
      };

  // Read the file
  size_t line_number = 0;
  std::string line, word;
  while (std::getline(file, line))
  {
    std::istringstream line_stream(line);
    line_stream >> word;

    auto& f = file;
    auto& ls = line_stream;
    auto& ln = line_number;

    // Get initial info, setup
    if (word == "NUM_GROUPS")
    {
      line_stream >> n_groups;
      reinit();
    }
    if (word == "NUM_MOMENTS")
    {
      assert(n_groups > 0);
      line_stream >> n_moments;
      transfer_matrices.resize(n_moments);
      for (auto& transfer_matrix : transfer_matrices)
        transfer_matrix.resize(n_groups, std::vector<double>(n_groups));
    }
    if (word == "NUM_PRECURSORS")
    {
      assert(n_groups > 0);
      line_stream >> n_precursors;
      precursor_lambda.assign(n_precursors, 0.0);
      precursor_yield.assign(n_precursors, 0.0);

      chi_delayed.resize(n_groups);
      for (unsigned int g = 0; g < n_groups; ++g)
        chi_delayed[g].resize(n_precursors);
    }
    if (word == "DENSITY")
      line_stream >> density;

    // Parse basic cross sections
    if (word == "SIGMA_T_BEGIN")
      read_1d_data("SIGMA_T", sigma_t, n_groups, f, ls, ln);
    if (word == "SIGMA_A_BEGIN")
      read_1d_data("SIGMA_A", sigma_a, n_groups, f, ls, ln);
    if (word == "SIGMA_F_BEGIN")
      read_1d_data("SIGMA_F", sigma_f, n_groups, f, ls, ln);

    if (word == "NU_SIGMA_F_BEGIN")
      read_1d_data("NU_SIGMA_F", nu_sigma_f, n_groups, f, ls, ln);
    if (word == "NU_PROMPT_SIGMA_F_BEGIN")
      read_1d_data("NU_PROMPT_SIGMA_F",
                   nu_prompt_sigma_f, n_groups, f, ls, ln);
    if (word == "NU_DELAYED_SIGMA_F_BEGIN")
      read_1d_data("NU_DELAYED_SIGMA_F",
                   nu_delayed_sigma_f, n_groups, f, ls, ln);

    if (word == "NU_BEGIN")
      read_1d_data("NU", nu, n_groups, f, ls, ln);
    if (word == "NU_PROMPT_BEGIN")
      read_1d_data("NU_PROMPT", nu_prompt, n_groups, f, ls, ln);
    if (word == "NU_DELAYED_BEGIN")
      read_1d_data("NU_DELAYED", nu_delayed, n_groups, f, ls, ln);
    if (word == "BETA_BEGIN")
      read_1d_data("BETA", beta, n_groups, f, ls, ln);

    if (word == "CHI_BEGIN")
      read_1d_data("CHI", chi, n_groups, f, ls, ln);
    if (word == "CHI_PROMPT_BEGIN")
      read_1d_data("CHI_PROMPT", chi_prompt, n_groups, f, ls, ln);

    if (word == "VELOCITY_BEGIN" and inv_velocity.front() == 0.0)
    {
      read_1d_data("VELOCITY", inv_velocity, n_groups, f, ls, ln);
      for (auto& v: inv_velocity) v = 1.0 / v;
    }
    if (word == "INV_VELOCITY_BEGIN")
      read_1d_data("INV_VELOCITY", inv_velocity, n_groups, f, ls, ln);

    if (word == "DIFFUSION_COEFF_BEGIN")
      read_1d_data("DIFFUSION_COEFF", diffusion_coeff, n_groups, f, ls, ln);
    if (word == "BUCKLING_BEGIN")
      read_1d_data("BUCKLING", buckling, n_groups, f, ls, ln);

    // Read transfer matrix
    if (word == "TRANSFER_MOMENTS_BEGIN")
      read_transfer_matrix("TRANSFER_MOMENTS", transfer_matrices, f, ls, ln);

    // Read delayed neutron data
    if (n_precursors > 0)
    {
      if (word == "PRECURSOR_LAMBDA_BEGIN")
        read_1d_data("PRECURSOR_LAMBDA",
                     precursor_lambda, n_precursors, f, ls, ln);
      if (word == "PRECURSOR_YIELD_BEGIN")
        read_1d_data("PRECURSOR_YIELD",
                     precursor_yield, n_precursors, f, ls, ln);
      if (word == "CHI_DELAYED_BEGIN")
        read_chi_delayed("CHI_DELAYED", chi_delayed, f, ls, ln);
    }

    word = "";
    ++line_number;
  }//while open
  file.close();

  compute_scattering_from_transfers();
  reconcile_cross_sections();
  reconcile_fission_properties();
  compute_macroscopic_cross_sections();

  if (verbose)
    std::cout << "Cross Sections Details:\n"
              << "\t# of Groups    : " << n_groups << "\n"
              << "\t# of Precursors: " << n_precursors << "\n"
              << "\tFissile?       : " << is_fissile << "\n";

  std::cout << "Done reading cross-sections file \"" << file_name << "\".\n";

}
