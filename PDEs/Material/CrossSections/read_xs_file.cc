#include "cross_sections.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <numeric>


void CrossSections::read_xs_file(const std::string& file_name)
{
  std::cout << "Reading cross-section file \"" << file_name << "\"\n";

  //========== Clear old cross sections
  reset();

  //========== Open the file
  std::ifstream file;
  file.open(file_name);
  if (!file.is_open()) {
    std::stringstream err;
    err << "CrossSections::" << __FUNCTION__ << ": "
        << "Failed to open the cross section file \""
        << file_name << "\".\n";
    throw std::runtime_error(err.str());
  }

  //========== Book-Keeping
  bool found_groups = false;
  bool found_velocity = false;
  bool found_inv_velocity = false;

  //========== Read the file
  unsigned int line_number = 0;
  std::string line, word;
  while (std::getline(file, line))
  {
    std::istringstream line_stream(line);
    line_stream >> word;

    //========== Shorthand
    auto& f = file;
    auto& ls = line_stream;
    auto& ln = line_number;

    //========== Get initial info, setup
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
      nu_sigma_f.assign(n_groups, 0.0);
      nu_prompt_sigma_f.assign(n_groups, 0.0);
      nu_delayed_sigma_f.assign(n_groups, 0.0);
      inv_velocity.assign(n_groups, 0.0);
      diffusion_coeff.assign(n_groups, 0.0);
      found_groups = true;
    }
    if (word == "SCATTERING_ORDER")
    {
      line_stream >> scattering_order;
      if (found_groups)
      {
        transfer_matrices.resize(scattering_order + 1);
        for (unsigned int m = 0; m < transfer_matrices.size(); ++m)
          transfer_matrices[m].resize(n_groups, std::vector<double>(n_groups));
      }
      else
      {
        std::stringstream err;
        err << "CrossSections::" << __FUNCTION__ << ": "
            << "NUM_GROUPS must be specified before SCATTERING_ORDER "
            << "in the xs file.";
        throw std::runtime_error(err.str());
      }
    }
    if (word == "NUM_PRECURSORS")
    {
      line_stream >> n_precursors;
      has_precursors = (n_precursors > 0);
      precursor_lambda.assign(n_precursors, 0.0);
      precursor_yield.assign(n_precursors, 0.0);

      if (found_groups)
      {
        chi_delayed.resize(n_groups);
        for (unsigned int g = 0; g < n_groups; ++g)
          chi_delayed[g].resize(n_precursors);
      }
      else
      {
        std::stringstream err;
        err << "CrossSections::" << __FUNCTION__ << ": "
            << "NUM_GROUPS must be specified before NUM_PRECURSORS "
            << "in the xs file.";
        throw std::runtime_error(err.str());
      }
    }

    // Parse basic cross sections
    if (word == "SIGMA_T_BEGIN")
      read_cross_section("SIGMA_T", sigma_t, f, ls, ln);
    if (word == "SIGMA_A_BEGIN")
      read_cross_section("SIGMA_A", sigma_a, f, ls, ln);
    if (word == "SIGMA_F_BEGIN")
      read_cross_section("SIGMA_F", sigma_f, f, ls, ln);

    if (word == "NU_BEGIN")
      read_cross_section("NU", nu, f, ls, ln);
    if (word == "NU_PROMPT_BEGIN")
      read_cross_section("NU_PROMPT", nu_prompt, f, ls, ln);
    if (word == "NU_DELAYED_BEGIN")
      read_cross_section("NU_DELAYED", nu_delayed, f, ls, ln);

    if (word == "CHI_BEGIN")
      read_cross_section("CHI", chi, f, ls, ln);
    if (word == "CHI_PROMPT_BEGIN")
      read_cross_section("CHI_PROMPT", chi_prompt, f, ls, ln);

    if (word == "VELOCITY_BEGIN" and not found_inv_velocity)
    {
      read_cross_section("VELOCITY", inv_velocity, f, ls, ln);
      found_velocity = true;
    }
    if (word == "INV_VELOCITY_BEGIN")
    {
      read_cross_section("INV_VELOCITY", inv_velocity, f, ls, ln);
      found_inv_velocity = true;
    }

    if (word == "DIFFUSION_COEFF_BEGIN")
      read_cross_section("DIFFUSION_COEFF", diffusion_coeff, f, ls, ln);

    // Read transfer matrix
    if (word == "TRANSFER_MOMENTS_BEGIN")
      read_transfer_matrices(
          "TRANSFER_MATRICES", transfer_matrices, f, ls, ln);

    // Read delayed neutron data
    if (has_precursors)
    {
      if (word == "PRECURSOR_LAMBDA_BEGIN")
        read_precursor_property("PRECURSOR_LAMBDA", precursor_lambda, f, ls, ln);
      if (word == "PRECURSOR_YIELD_BEGIN")
        read_precursor_property("PRECURSOR_YIELD", precursor_yield, f, ls, ln);
      if (word == "CHI_DELAYED_BEGIN")
        read_delayed_spectra("CHI_DELAYED", chi_delayed, f, ls, ln);
    }

    word = ""; ++line_number;
  }//while open
  file.close();

  // Convert velocity to inverse velocity
  if (found_velocity)
    for (auto& v : inv_velocity)
      v = 1.0 / v;

  // Compute the scattering xs from the transfer matrix
  compute_scattering_from_transfers();

  // Enforce realistic total, abosrption, and scattering xs
  reconcile_cross_sections();

  // Reconcile the fission properties
  reconcile_fission_properties();

  // Compute macroscopic cross sections
  compute_macroscopic_cross_sections();

  std::cout << "Finished reading cross-section file \"" << file_name << "\"\n";
}
