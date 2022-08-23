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

  // Book-Keeping
  bool found_groups = false;
  bool found_inv_velocity = false;

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
      sigma_t.assign(n_groups, 0.0);
      sigma_a.assign(n_groups, 0.0);
      sigma_s.assign(n_groups, 0.0);
      sigma_r.assign(n_groups, 0.0);
      sigma_f.assign(n_groups, 0.0);
      chi.assign(n_groups, 0.0);
      chi_prompt.assign(n_groups, 0.0);
      nu.assign(n_groups, 0.0);
      nu_prompt.assign(n_groups, 0.0);
      nu_delayed.assign(n_groups, 0.0);
      beta.assign(n_groups, 0.0);
      nu_sigma_f.assign(n_groups, 0.0);
      nu_prompt_sigma_f.assign(n_groups, 0.0);
      nu_delayed_sigma_f.assign(n_groups, 0.0);
      inv_velocity.assign(n_groups, 0.0);
      diffusion_coeff.assign(n_groups, 0.0);
      buckling.assign(n_groups, 0.0);
      found_groups = true;
    }
    if (word == "SCATTERING_ORDER")
    {
      assert(found_groups);

      line_stream >> scattering_order;
      transfer_matrices.resize(scattering_order + 1);
      for (auto& transfer_matrix: transfer_matrices)
        transfer_matrix.resize(n_groups, std::vector<double>(n_groups));
    }
    if (word == "NUM_PRECURSORS")
    {
      assert(found_groups);
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
      read_cross_section("SIGMA_T", sigma_t, f, ls, ln);
    if (word == "SIGMA_A_BEGIN")
      read_cross_section("SIGMA_A", sigma_a, f, ls, ln);
    if (word == "SIGMA_F_BEGIN")
    {
      read_cross_section("SIGMA_F", sigma_f, f, ls, ln);
      is_fissile = (is_fissile) ||
                   std::accumulate(sigma_f.begin(), sigma_f.end(), 0.0) > 1.0e-12;
    }
    if (word == "NU_SIGMA_F_BEGIN")
    {
      read_cross_section("NU_SIGMA_F", nu_sigma_f, f, ls, ln);
      is_fissile = (is_fissile) ||
                   std::accumulate(nu_sigma_f.begin(), nu_sigma_f.end(), 0.0) > 1.0e-12;
    }
    if (word == "NU_PROMPT_SIGMA_F_BEGIN")
    {
      read_cross_section("NU_PROMPT_SIGMA_F", nu_prompt_sigma_f, f, ls, ln);
      is_fissile = (is_fissile) ||
                   std::accumulate(nu_prompt_sigma_f.begin(),
                                   nu_prompt_sigma_f.end(), 0.0) > 1.0e-12;
    }
    if (word == "NU_DELAYED_SIGMA_F_BEGIN")
    {
      read_cross_section("NU_DELAYED_SIGMA_F", nu_delayed_sigma_f, f, ls, ln);
      is_fissile = (is_fissile) ||
                   std::accumulate(nu_delayed_sigma_f.begin(),
                                   nu_delayed_sigma_f.end(), 0.0) > 1.0e-12;
    }

    if (word == "NU_BEGIN")
      read_cross_section("NU", nu, f, ls, ln);
    if (word == "NU_PROMPT_BEGIN")
      read_cross_section("NU_PROMPT", nu_prompt, f, ls, ln);
    if (word == "NU_DELAYED_BEGIN")
      read_cross_section("NU_DELAYED", nu_delayed, f, ls, ln);
    if (word == "BETA_BEGIN")
      read_cross_section("BETA", beta, f, ls, ln);

    if (word == "CHI_BEGIN")
      read_cross_section("CHI", chi, f, ls, ln);
    if (word == "CHI_PROMPT_BEGIN")
      read_cross_section("CHI_PROMPT", chi_prompt, f, ls, ln);

    if (word == "VELOCITY_BEGIN" and !found_inv_velocity)
    {
      read_cross_section("VELOCITY", inv_velocity, f, ls, ln);
      for (auto& v: inv_velocity) v = 1.0 / v;
    }
    if (word == "INV_VELOCITY_BEGIN")
    {
      read_cross_section("INV_VELOCITY", inv_velocity, f, ls, ln);
      found_inv_velocity = true;
    }

    if (word == "DIFFUSION_COEFF_BEGIN")
      read_cross_section("DIFFUSION_COEFF", diffusion_coeff, f, ls, ln);
    if (word == "BUCKLING_BEGIN")
      read_cross_section("BUCKLING", buckling, f, ls, ln);

    // Read transfer matrix
    if (word == "TRANSFER_MOMENTS_BEGIN")
      read_transfer_matrices(
          "TRANSFER_MOMENTS", transfer_matrices, f, ls, ln);

    // Read delayed neutron data
    if (n_precursors > 0)
    {
      if (word == "PRECURSOR_LAMBDA_BEGIN")
        read_precursor_property("PRECURSOR_LAMBDA", precursor_lambda, f, ls, ln);
      if (word == "PRECURSOR_YIELD_BEGIN")
        read_precursor_property("PRECURSOR_YIELD", precursor_yield, f, ls, ln);
      if (word == "CHI_DELAYED_BEGIN")
        read_delayed_spectra("CHI_DELAYED", chi_delayed, f, ls, ln);
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
