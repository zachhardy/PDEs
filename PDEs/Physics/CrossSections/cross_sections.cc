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
  for (size_t gp = 0; gp < n_groups; ++gp)
  {
    for (size_t g = 0; g < n_groups; ++g)
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
    for (size_t g = 0; g < n_groups; ++g)
    {
      assert(sigma_t[g] >= 0.0);
      assert(sigma_t[g] >= sigma_s[g]);
      sigma_a[g] = sigma_t[g] - sigma_s[g];
    }
  else
    for (size_t g = 0; g < n_groups; ++g)
    {
      assert(sigma_a[g] >= 0.0);
      sigma_t[g] = sigma_a[g] + sigma_s[g];
  }

  // Compute the removal cross sections
  for (size_t g = 0; g < n_groups; ++g)
    sigma_r[g] = sigma_t[g] - transfer_matrices[0][g][g];
}


//######################################################################


void
Physics::CrossSections::reconcile_fission_properties()
{
  // Check whether the material is fissile
  double sum_sigma_f = std::accumulate(sigma_f.begin(), sigma_f.end(), 0.0);
  is_fissile = (sum_sigma_f > 1.0e-12);

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
    // Check for negative cross sections
    for (size_t g = 0; g < n_groups; ++g)
    assert(sigma_f[g] >= 0.0);

    // Determine which terms are present
    std::pair<bool, bool> has_total(false, false);
    std::pair<bool, bool> has_prompt(false, false);
    std::pair<bool, bool> has_delayed(false, false);
    for (size_t g = 0; g < n_groups; ++g)
    {
      // Total fission properties
      if (not has_total.first and nu[g] > 0.0)
        has_total.first = true;
      if (not has_total.second and chi[g] > 0.0)
        has_total.second = true;

      // Prompt fission properties
      if (not has_prompt.first and nu_prompt[g] > 0.0)
        has_prompt.first = true;
      if (not has_prompt.second and chi_prompt[g] > 0.0)
        has_prompt.second = true;

      // Delayed fission properties
      if (not has_delayed.first and nu_delayed[g] > 0.0)
        has_delayed.first = true;
    }

    // Check for prompt/delayed quantities
    if (n_precursors > 0)
    {
      assert(has_prompt.first && has_prompt.second && has_delayed.first);

      // Ensure positive terms
      for (size_t g = 0; g < n_groups; ++g)
      {
        assert(nu_prompt[g] >= 0.0 && nu_delayed[g] >= 0.0);
        assert(chi_prompt[g] >= 0.0);
      }

      // Check precursor properties
      for (size_t j = 0; j < n_precursors; ++j)
      {
        assert(precursor_lambda[j] > 0.0);
        assert(precursor_yield[j] > 0.0);
      }

      // Check delayed spectra
      for (size_t j = 0; j < n_precursors; ++j)
      {
        // Ensure positivity
        for (size_t g = 0; g < n_groups; ++g)
          assert(chi_delayed[g][j] >= 0.0);

        // Compute spectra sum
        double sum = 0.0;
        for (size_t g = 0; g < n_groups; ++g)
          sum += chi_delayed[g][j];
        assert(sum > 0.0);
      }
      has_delayed.second = true;

      // Normalize prompt spectra
      double prompt_sum = std::accumulate(chi_prompt.begin(),
                                          chi_prompt.end(), 0.0);

      if (std::abs(prompt_sum - 1.0) > 1.0e-12)
      {
        std::stringstream warn;
        warn << "!!! Warning !!! CrossSections::" << __FUNCTION__ << ": "
             << "Normalizing prompt fission spectrum to unity.";
        std::cout << warn.str() << std::endl;

        for (auto& v : chi_prompt)
          v /= prompt_sum;
      }

      // Normalize delayed spectra
      for (size_t j = 0; j < n_precursors; ++j)
      {
        double sum = 0.0;
        for (size_t g = 0; g < n_groups; ++g)
          sum += chi_delayed[g][j];

        if (std::abs(sum - 1.0) > 1.0e-12)
        {
          std::stringstream warn;
          warn << "!!! Warning !!! CrossSections::" << __FUNCTION__ << ": "
               << "Normalizing delayed emmission spectrum to unity "
               << "for precursor species " << j << ".";
          std::cout << warn.str() << std::endl;

          for (size_t g = 0; g < n_groups; ++g)
            chi_delayed[g][j] /= sum;
        }
      }

      // Normalize precursor yields
      double yield_sum = std::accumulate(precursor_yield.begin(),
                                         precursor_yield.end(), 0.0);
      if (std::abs(yield_sum - 1.0) > 1.0e-12)
      {
        std::stringstream warn;
        warn << "!!! Warning !!! CrossSections::" << __FUNCTION__ << ": "
             << "Normalizing precursor yields to unity.";
        std::cout << warn.str() << std::endl;

        for (auto& v : precursor_yield)
          v /= yield_sum;
      }

      // Compute total quantities
      for (size_t g = 0; g < n_groups; ++g)
      {
        nu[g] = nu_prompt[g] + nu_delayed[g];

        // Compute the \beta-weighted total fission spectra. In this context
        // \beta is the fraction of fissions that are delayed. The prompt
        // fission spectra is weighted by `1 - \beta` and the delayed emmission
        // spectra by \beta \gamma_j, where \gamma_j is the yield of precursor
        // species j.
        double beta = nu_delayed[g]/nu[g];
        chi[g] = (1.0 - beta)*chi_prompt[g];
        for (size_t j = 0; j < n_precursors; ++j)
          chi[g] += beta*precursor_yield[j]*chi_delayed[g][j];
      }
      has_total.first = true;
      has_total.second = true;
    }//if has_precursors
    else
    {
      assert(has_total.first && has_total.second);
      assert(has_total.first && has_total.second);

      // Normalize fission spectrum
      double sum = std::accumulate(chi.begin(), chi.end(), 0.0);
      if (std::abs(sum - 1.0) > 1.0e-12)
      {
        std::stringstream warn;
        warn << "!!! Warning !!! CrossSections::" << __FUNCTION__ << ": "
             << "Normalizing total fission spectrum to unity.";
        std::cout << warn.str() << std::endl;

        for (auto& v : chi)
          v /= sum;
      }
    }

    // Compute nu_sigma_f terms
    for (size_t g = 0; g < n_groups; ++g)
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
  for (size_t g = 0; g < n_groups; ++g)
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

  for (size_t m = 0; m <= scattering_order; ++m)
    for (size_t g = 0; g < n_groups; ++g)
      for (size_t gp = 0; gp < n_groups; ++gp)
        transfer_matrices[m][g][gp] *= density;

  // Compute diffusion coefficient if unspecified
  double sum = std::accumulate(diffusion_coeff.begin(),
                               diffusion_coeff.end(), 0.0);
  if (sum < 1.0e-12)
  {
    for (size_t g = 0; g < n_groups; ++g)
      diffusion_coeff[g] = 1.0/(3.0*sigma_t[g]);
  }
}
