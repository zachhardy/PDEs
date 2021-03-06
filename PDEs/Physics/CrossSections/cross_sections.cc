#include "cross_sections.h"
#include "macros.h"

#include <iostream>
#include <sstream>
#include <cmath>
#include <numeric>


Physics::CrossSections::CrossSections() :
  MaterialProperty(MaterialPropertyType::CROSS_SECTIONS)
{}


Physics::CrossSections::CrossSections(const std::string property_name) :
  MaterialProperty(MaterialPropertyType::CROSS_SECTIONS, property_name)
{}

//######################################################################


void
Physics::CrossSections::reset()
{
  n_groups = 0;
  n_precursors = 0;
  scattering_order = 0;

  is_fissile = false;
  density = 1.0;

  sigma_t.clear();
  sigma_a.clear();
  sigma_s.clear();
  sigma_r.clear();
  sigma_f.clear();

  chi.clear();
  chi_prompt.clear();
  chi_delayed.clear();

  nu.clear();
  nu_prompt.clear();
  nu_delayed.clear();

  nu_sigma_f.clear();
  nu_prompt_sigma_f.clear();
  nu_delayed_sigma_f.clear();

  precursor_lambda.clear();
  precursor_yield.clear();

  inv_velocity.clear();
  diffusion_coeff.clear();

  transfer_matrices.clear();
}


//######################################################################


void
Physics::CrossSections::compute_scattering_from_transfers()
{
  sigma_s.assign(n_groups, 0.0);
  if (transfer_matrices.empty())
    return;

  /* A group's scattering cross section is defined as the sum of the
   * transfer cross sections from a fixed group to all  other groups. In the
   * transfer matrices, the rows contain the destination groups and columns
   * the origin group. Due to this, computing sigma_s necessitates a column-
   * wise sum. */
  for (unsigned int gp = 0; gp < n_groups; ++gp)
  {
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      assert(transfer_matrices[0][g][gp] >= 0.0);
      sigma_s[gp] += transfer_matrices[0][g][gp];
    }
  }
}


//######################################################################


void
Physics::CrossSections::reconcile_cross_sections()
{

  // Determine whether sigma_a was specified
  double sum = std::accumulate(sigma_a.begin(), sigma_a.end(), 0.0);
  bool specified_sigma_a = (sum > 1.0e-12);

  /* Compute absorption xs from transfer matrix, if not specified. Otherwise,
   * recompute the total cross-section from the specified scattering and
   * absorption. */
  if (not specified_sigma_a)
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      assert(sigma_t[g] >= 0.0);
      assert(sigma_t[g] >= sigma_s[g]);
      sigma_a[g] = sigma_t[g] - sigma_s[g];
    }
  else
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      assert(sigma_a[g] >= 0.0);
      sigma_t[g] = sigma_a[g] + sigma_s[g];
  }

  // Compute the removal cross sections
  for (unsigned int g = 0; g < n_groups; ++g)
    sigma_r[g] = sigma_t[g] - transfer_matrices[0][g][g];
}


//######################################################################


void
Physics::CrossSections::reconcile_fission_properties()
{
  auto check_xs = [](std::vector<double> x)
  {
    for (const auto& v : x)
      if (v != 0.0)
        return true;
    return false;
  };

  auto check_matrix = [](std::vector<std::vector<double>> A)
  {
    for (const auto& a : A)
      for (const auto& v : a)
        if (v != 0.0)
          return true;
    return false;
  };

  // Clear precursors if not fissile
  if (not is_fissile and n_precursors > 0)
  {
    n_precursors = 0;
    precursor_lambda.clear();
    precursor_yield.clear();
    chi_delayed.clear();

    std::cout << "!!! WARNING !!! " << "CrossSections::" << __FUNCTION__
              << ": Precursors found in a non-fissile material.\n"
              << "Clearing the precursor properties for consistency.";
  }

  // Check fission properties
  if (is_fissile)
  {
    // Check which quantities were specified
    const bool has_sigf = check_xs(sigma_f);
    const bool has_nusigf = check_xs(nu_sigma_f);
    const bool has_nupsigf = check_xs(nu_prompt_sigma_f);
    const bool has_nudsigf = check_xs(nu_delayed_sigma_f);

    const bool has_nu = check_xs(nu);
    const bool has_nup = check_xs(nu_prompt);
    const bool has_nud = check_xs(nu_delayed);
    const bool has_beta = check_xs(beta);

    const bool has_chi = check_xs(chi);
    const bool has_chip = check_xs(chi_prompt);
    const bool has_chid = check_matrix(chi_delayed);

    /* Check the specified properties from most general to least. Preference
     * is given when the explicit properties are provided. For example, if
     * prompt and delayed quantities are explicitly specified, the values are
     * kept and others derived from them. If explicit prompt and delayed
     * quantities are not specified, alternate specification forms are checked.
     * At present, the only alternate way to specify prompt and delayed
     * quantities is to use the ``delayed fraction'' method where nu_sigma_f,
     * nu, and beta are specified group-wise. From these, all prompt and
     * delayed quantities can be derived. If neither of these conditions are
     * met, it is assumed that only total fission is desired and the prompt and
     * delayed properties are set accordingly. */

    // Begin with prompt/delayed quantities
    if (n_precursors > 0)
    {
      // Full specification
      if (has_sigf and has_nup and has_nud)
      {
        assert(std::all_of(sigma_f.begin(), sigma_f.end(),
                           [](double x) { return x >= 0.0; }));
        assert(std::all_of(nu_prompt.begin(), nu_prompt.end(),
                           [](double x) { return x > 0.0; }));
        assert(std::all_of(nu_delayed.begin(), nu_delayed.end(),
                           [](double x) { return x >= 0.0; }));

        for (unsigned int g = 0; g < n_groups; ++g)
          beta[g] = nu_delayed[g] / nu[g];
      }

      // Delayed-fraction specification
      else if (has_nusigf and has_nu and has_beta)
      {
        assert(std::all_of(nu_sigma_f.begin(), nu_sigma_f.end(),
                           [](double x) { return x >= 0.0; }));
        assert(std::all_of(nu.begin(), nu.end(),
                           [](double x) { return x >= 1.0; }));
        assert(std::all_of(beta.begin(), beta.end(),
                           [](double x) { return x >= 0.0; }));

        for (unsigned int g = 0; g < n_groups; ++g)
        {
          sigma_f[g] = nu_sigma_f[g] / nu[g];
          nu_prompt[g] = (1.0 - beta[g]) * nu[g];
          nu_delayed[g] = beta[g] * nu[g];
        }
      }

      // Check the spectra
      assert(has_chip && has_chid);
      assert(std::all_of(chi_prompt.begin(), chi_prompt.end(),
                         [](double x) { return x >= 0.0; }));
      assert(std::all_of(chi_delayed.begin(), chi_delayed.end(),
                         [](std::vector<double> x)
                         { return std::all_of(x.begin(), x.end(),
                                              [](double y)
                                              { return y >= 0.0; }); }));

      // Check precursor properties
      assert(std::all_of(precursor_yield.begin(), precursor_yield.end(),
                         [](double x) { return x >= 0.0; }));
      assert(std::all_of(precursor_lambda.begin(), precursor_lambda.end(),
                         [](double x) { return x > 0.0; }));

      // Normalize the prompt spectra
      double prompt_sum =
          std::accumulate(chi_prompt.begin(), chi_prompt.end(), 0.0);

      if (std::abs(prompt_sum - 1.0) > 1.0e-12)
      {
        std::stringstream warn;
        warn << "**!!!** Warning **!!!** CrossSections::"
             << __FUNCTION__ << ": "
             << "Normalizing prompt fission spectrum to unity.";
        std::cout << warn.str() << std::endl;

        for (unsigned int g = 0; g < n_groups; ++g)
          chi_prompt[g] /= prompt_sum;
      }

      // Normalize the delayed spectra
      for (unsigned int j = 0; j < n_precursors; ++j)
      {
        double delayed_sum_j = 0.0;
        for (unsigned int g = 0; g < n_groups; ++g)
          delayed_sum_j += chi_delayed[g][j];

        if (std::abs(delayed_sum_j - 1.0) > 1.0e-12)
        {
          std::stringstream warn;
          warn << "**!!!** Warning **!!!** CrossSections::"
               << __FUNCTION__ << ": "
               << "Normalizing delayed emission spectrum to unity "
               << "for precursor species " << j << ".";
          std::cout << warn.str() << std::endl;
        }

        for (unsigned int g = 0; g < n_groups; ++g)
          chi_delayed[g][j] /= delayed_sum_j;
      }

      // Normalize the precursor yields
      double yield_sum =
          std::accumulate(precursor_yield.begin(), precursor_yield.end(), 0.0);

      if (std::abs(yield_sum - 1.0) > 1.0e-12)
      {
        std::stringstream warn;
        warn << "**!!!** Warning **!!!** CrossSections::"
             << __FUNCTION__ << ": "
             << "Normalizing precursor yields to unity.";
        std::cout << warn.str() << std::endl;

        for (unsigned int j = 0; j < n_precursors; ++j)
          precursor_yield[j] /= yield_sum;
      }

      // Compute the total fission quantities from prompt and delayed
      for (unsigned int g = 0; g < n_groups; ++g)
      {
        nu[g] = nu_prompt[g] + nu_delayed[g];

        // Compute the beta-weighted total fission spectra
        chi[g] = (1.0 - beta[g])*chi_prompt[g];
        for (unsigned int j = 0; j < n_precursors; ++j)
          chi[g] += beta[g]*precursor_yield[j]*chi_delayed[g][j];
      }

      // Check values
      assert(std::all_of(nu.begin(), nu.end(),
                         [](double x) { return x >= 1.0; }));
      assert(std::all_of(chi.begin(), chi.end(),
                         [](double x) { return x >= 0.0; }));
      double chi_sum = std::accumulate(chi.begin(), chi.end(), 0.0);
      assert(std::abs(chi_sum - 1.0) < 1.0e-12);
    }//if precursors

    // Total fission only
    else
    {
      assert(has_chi && has_nu && (has_nusigf || has_sigf));
      assert(std::all_of(nu.begin(), nu.end(),
                         [](double x) { return x >= 1.0; }));
      assert(std::all_of(chi.begin(), chi.end(),
                         [](double x) { return x >= 0.0; }));

      if (has_nusigf and has_nu and not has_sigf)
        for (unsigned int g = 0; g < n_groups; ++g)
          sigma_f[g] = nu_sigma_f[g] / nu[g];

      assert(std::all_of(sigma_f.begin(), sigma_f.end(),
                         [](double x) { return x >= 0.0; }));

      // Normalize the fission spectrum
      double sum = std::accumulate(chi.begin(), chi.end(), 0.0);
      if (std::abs(sum - 1.0) > 1.0e-12)
      {
        std::stringstream warn;
        warn << "**!!!** Warning **!!!** CrossSections::"
             << __FUNCTION__ << ": "
             << "Normalizing total fission spectrum to unity.";
        std::cout << warn.str() << std::endl;

        for (unsigned int g = 0; g < n_groups; ++g)
          chi[g] /= sum;
      }
    }

    // Compute nu_sigma_f terms
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      nu_sigma_f[g] = nu[g]*sigma_f[g];
      nu_prompt_sigma_f[g] = nu_prompt[g]*sigma_f[g];
      nu_delayed_sigma_f[g] = nu_delayed[g]*sigma_f[g];
    }
  }//if fissile
}


//######################################################################


void
Physics::CrossSections::compute_macroscopic_cross_sections()
{
  for (unsigned int g = 0; g < n_groups; ++g)
  {
    sigma_t[g] *= density;
    sigma_a[g] *= density;
    sigma_s[g] *= density;
    sigma_r[g] *= density;

    sigma_f[g] *= density;
    nu_sigma_f[g] *= density;
    nu_prompt_sigma_f[g] *= density;
    nu_delayed_sigma_f[g] *= density;
  }

  for (unsigned int m = 0; m <= scattering_order; ++m)
    for (unsigned int g = 0; g < n_groups; ++g)
      for (unsigned int gp = 0; gp < n_groups; ++gp)
        transfer_matrices[m][g][gp] *= density;

  // Compute diffusion coefficient if unspecified
  double sum = std::accumulate(diffusion_coeff.begin(),
                               diffusion_coeff.end(), 0.0);
  if (sum < 1.0e-12)
  {
    for (unsigned int g = 0; g < n_groups; ++g)
      diffusion_coeff[g] = 1.0/(3.0*sigma_t[g]);
  }
}
