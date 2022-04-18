#include "cross_sections.h"

#include <iostream>
#include <sstream>
#include <numeric>


//######################################################################
/**
 * Compute the group-wise scattering cross sections from the zeroth transfer
 * matrices. This is defined as the sum of all transfers from a fixed group to
 * any other. Mathematically, this is given by
 * \f[ \sigma_{s,g} = \sum_{g^\prime} \sigma_{0, g \rightarrow g^\prime} ,\f]
 * which is obtained via column-wise sums.
 */
void CrossSections::compute_scattering_from_transfers()
{
  sigma_s.assign(n_groups, 0.0);
  if (transfer_matrices.empty())
    return;

  // A group's scattering cross section is defined as the sum of the
  // transfer cross sections from a fixed group to all  other groups. In the
  // transfer matrices, the rows contain the destination groups and columns
  // the origin group. Due to this, computing sigma_s necessitates a column-
  // wise sum.
  for (unsigned int gp = 0; gp < n_groups; ++gp)
  {
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      if (transfer_matrices[0][g][gp] < 0.0)
      {
        std::stringstream err;
        err << "CrossSections::" << __FUNCTION__ << ": "
            << "Negative transfer cross section encountered for "
            << "starting group " << gp << " and destination group"
            << g << "." ;
        throw std::runtime_error(err.str());
      }
      sigma_s[gp] += transfer_matrices[0][g][gp];
    }
  }
}


//######################################################################
/**
 * If the absorption cross section was not specified, then compute it
 * via \f$ \sigma_a = \sigma_t - \sigma_s , \f$ where \f$ \sigma_s \f$ is
 * obtained from the transfer matrix. Otherwise, modify \f$ \sigma_t \f$ using
 * the provided absorption cross section and scattering cross section obtained
 * from the transfer matrix. If both \f$ \sigma_t \f$ and \f$ \sigma_a \f$ are
 * provided and they do not agree with the transfer matrix, the \f$ \sigma_a \f$
 * values are taken as true.
 */
void CrossSections::reconcile_cross_sections()
{

  // Determine whether sigma_a was specified
  double sum = std::accumulate(sigma_a.begin(), sigma_a.end(), 0.0);
  bool specified_sigma_a = (sum > 1.0e-12);

  // Compute aborption xs from transfer matrix, if not specified
  if (not specified_sigma_a)
  {
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      // Ensure positivity
      if (sigma_t[g] < 0.0)
      {
        std::stringstream err;
        err << "CrossSections::" << __FUNCTION__ << ": "
            << "Negative total cross section encountered in "
            << "group " << g << "." ;
        throw std::runtime_error(err.str());
      }

      // Ensure physically realistic setup
      if (sigma_t[g] < sigma_s[g])
      {
        std::stringstream err;
        err << "CrossSections::" << __FUNCTION__ << ": "
            << "Scattering cross section is larger than total for "
            << "group " << g << "." ;
        throw std::runtime_error(err.str());
      }

      sigma_a[g] = sigma_t[g] - sigma_s[g];
    }
  }

  // Otherwise, recompute the total xs from absorption and scattering
  else
  {
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      // Ensure positivity
      if (sigma_a[g] < 0.0)
      {
        std::stringstream err;
        err << "CrossSections::" << __FUNCTION__ << ": "
            << "Negative absorption cross section encountered in "
            << "group " << g << "." ;
        throw std::runtime_error(err.str());
      }

      sigma_t[g] = sigma_a[g] + sigma_s[g];
    }
  }

  // Compute the removal cross sections
  for (unsigned int g = 0; g < n_groups; ++g)
    sigma_r[g] = sigma_t[g] - transfer_matrices[0][g][g];
}


//######################################################################
/**
 * This routine does a number of things.
 * 1. If non-fissile but delayed neutron precursor properties were specified,
 *    these are cleared.
 * 2. If fissile and delayed neutron precursor properties were specified, checks
 *    for prompt and delayed \f$ \nu \f$ and \f$ \chi \f$ are performed and
 *    fission/emmission spectra are normalized to unity. Additionally, precursor
 *    yields \f$ \gamma \f$ are normalized to unity. Lastly, the total
 *    \f$ \nu \f$ and \f$ \chi \f$ quantities are computed from their prompt
 *    and delayed counterparts.
 * 3. If fissile and no delayed neutron precursor properties were specified,
 *    checks for total \f$ \nu \f$ and \f$ \chi \f$ are performed and the
 *    fission spectrum is normalized to unity.
 */
void CrossSections::reconcile_fission_properties()
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
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      if (sigma_f[g] < 0.0)
      {
        std::stringstream err;
        err << "CrossSections::" << __FUNCTION__ << ": "
            << "Negative fisson cross section encountered in "
            << "group " << g << "." ;
        throw std::runtime_error(err.str());
      }
    }

    // Determine which terms are present
    std::pair<bool, bool> has_total(false, false);
    std::pair<bool, bool> has_prompt(false, false);
    std::pair<bool, bool> has_delayed(false, false);
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      if (not has_total.first and nu[g] > 0.0)
        has_total.first = true;
      if (not has_total.second and chi[g] > 0.0)
        has_total.second = true;

      if (not has_prompt.first and nu_prompt[g] > 0.0)
        has_prompt.first = true;
      if (not has_prompt.second and chi_prompt[g] > 0.0)
        has_prompt.second = true;

      if (not has_delayed.first and nu_delayed[g] > 0.0)
        has_delayed.first = true;
    }

    // Check for prompt/delayed problems
    if (has_precursors)
    {
      // Ensure prompt/delayed quantities are provided
      if (not has_prompt.first or
          not has_prompt.second or
          not has_delayed.first)
      {
        std::stringstream err;
        err << "CrossSections::" << __FUNCTION__ << ": "
            << "Prompt and delayed fission terms must be supplied for  "
            << "problems with delayed neutron precursors.";
        throw std::runtime_error(err.str());
      }

      // Ensure positivity
      for (unsigned int g = 0; g < n_groups; ++g)
      {
        if (nu_prompt[g] < 0.0 or nu_delayed[g] < 0.0)
        {
          std::stringstream err;
          err << "CrossSections::" << __FUNCTION__ << ": "
              << "Negative nu prompt or nu delayed value encountered in "
              << "group " << g << "." ;
          throw std::runtime_error(err.str());
        }
        if (chi_prompt[g] < 0.0)
        {
          std::stringstream err;
          err << "CrossSections::" << __FUNCTION__ << ": "
              << "Negative prompt fission spectrum value encountered in "
              << "group " << g << "." ;
          throw std::runtime_error(err.str());
        }
      }

      // Check precursor properties
      for (unsigned int j = 0; j < n_precursors; ++j)
      {
        if (precursor_lambda[j] < 1.0e-12)
        {
          std::stringstream err;
          err << "CrossSections::" << __FUNCTION__ << ": "
              << "Decay constants must be non-zero.\n"
              << "Issue with precursor species " << j << ".";
          throw std::runtime_error(err.str());
        }
        if (precursor_yield[j] < 1.0e-12)
        {
          std::stringstream err;
          err << "CrossSections::" << __FUNCTION__  << ": "
              << "Yields must be non-zero.\n"
              << "Issue with precursor species " << j << ".";
          throw std::runtime_error(err.str());
        }
      }

      // Check delayed spectra
      for (unsigned int j = 0; j < n_precursors; ++j)
      {
        // Ensure positivity
        for (unsigned int g = 0; g < n_groups; ++g)
        {
          if (chi_delayed[g][j] < 0.0)
          {
            std::stringstream err;
            err << "CrossSections::" << __FUNCTION__ << ": "
                << "Negative delayed emmission spectrum value encountered in "
                << "group " << g << " for precursor species " << j << "." ;
            throw std::runtime_error(err.str());
          }
        }

        // Compute spectra sum
        double sum = 0.0;
        for (unsigned int g = 0; g < n_groups; ++g)
          sum += chi_delayed[g][j];

        if (sum < 1.0e-12)
        {
          std::stringstream err;
          err << "CrossSections::" << __FUNCTION__ << ": "
              << "Emmission spectra must be non-zero.\n"
              << "Issue encountered with precursor species " << j << ".";
          throw std::runtime_error(err.str());
        }
      }
      has_delayed.second = true;

      // Normalize prompt spectra
      double prompt_sum = std::accumulate(chi_prompt.begin(),
                                          chi_prompt.end(), 0.0);
      if (std::abs(prompt_sum - 1.0) > 1.0e-12)
      {
        std::stringstream warn;
        warn  << "!!! Warning !!! CrossSections::" << __FUNCTION__ << ": "
              << "Normalizing prompt fission spectrum to unity.";
        std::cout << warn.str() << std::endl;

        for (auto& v : chi_prompt) v /= prompt_sum;
      }

      // Normalize delayed spectra
      for (unsigned int j = 0; j < n_precursors; ++j)
      {
        double sum = 0.0;
        for (unsigned int g = 0; g < n_groups; ++g)
          sum += chi_delayed[g][j];

        if (std::abs(sum - 1.0) > 1.0e-12)
        {
          std::stringstream warn;
          warn  << "!!! Warning !!! CrossSections::" << __FUNCTION__ << ": "
                << "Normalizing delayed emmission spectrum to unity "
                << "for precursor species " << j << ".";
          std::cout << warn.str() << std::endl;

          for (unsigned int g = 0; g < n_groups; ++g)
            chi_delayed[g][j] /= sum;
        }
      }

      // Normalize precursor yields
      double yield_sum = std::accumulate(precursor_yield.begin(),
                                         precursor_yield.end(), 0.0);
      if (std::abs(yield_sum - 1.0) > 1.0e-12)
      {
        std::stringstream warn;
        warn  << "!!! Warning !!! CrossSections::" << __FUNCTION__ << ": "
              << "Normalizing precursor yields to unity.";
        std::cout << warn.str() << std::endl;

        for (auto& v : precursor_yield) v /= yield_sum;
      }

      // Compute total quantities
      for (unsigned int g = 0; g < n_groups; ++g)
      {
        nu[g] = nu_prompt[g] + nu_delayed[g];

        // Compute the \beta-weighted total fission spectra. In this context
        // \beta is the fraction of fissions that are delayed. The prompt
        // fission spectra is weighted by `1 - \beta` and the delayed emmission
        // spectra by \beta \gamma_j, where \gamma_j is the yield of precursor
        // species j.
        double beta = nu_delayed[g] / nu[g];
        chi[g] = (1.0 - beta) * chi_prompt[g];
        for (unsigned int j = 0; j < n_precursors; ++j)
          chi[g] = beta * precursor_yield[j] * chi_delayed[g][j];
      }
      has_total.first = true;
      has_total.second= true;
    }//if has_precursors

    // Handle problems w/out precursors
    else
    {
      // Ensure the total quantities are provided
      if (not has_total.first or not has_total.second)
      {
        std::stringstream err;
        err << "CrossSections::" << __FUNCTION__ << ": "
            << "Total fission terms were not found nor computed from "
            << "prompt and delayed fission terms.";
        throw std::runtime_error(err.str());
      }

      // Normalize fission spectrum
      double sum = std::accumulate(chi.begin(), chi.end(), 0.0);
      if (std::abs(sum - 1.0) > 1.0e-12)
      {
        std::stringstream warn;
        warn  << "!!! Warning !!! CrossSections::" << __FUNCTION__ << ": "
              << "Normalizing total fission spectrum to unity.";
        std::cout << warn.str() << std::endl;

        for (auto& v : chi) v /= sum;
      }
    }

    // Compute nu_sigma_f terms
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      nu_sigma_f[g] = nu[g] * sigma_f[g];
      nu_prompt_sigma_f[g] = nu_prompt[g] * sigma_f[g];
      nu_delayed_sigma_f[g] = nu_delayed[g] * sigma_f[g];
    }
  }//if fissile
}


//######################################################################
/**
 * Compute the macroscopic cross sections via
 * \f[ \Sigma_x = \rho \sigma_x .\f]
 * If the \p diffusion_coeff was unspecified, it is computed via its standard
 * defintion.
 */
void CrossSections::compute_macroscopic_cross_sections()
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
      diffusion_coeff[g] = 1.0 / (3.0 * sigma_t[g]);
  }
}
