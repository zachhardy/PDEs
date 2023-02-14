#include "cross_sections.h"
#include "macros.h"

#include <iostream>
#include <sstream>
#include <cmath>
#include <numeric>

using namespace std;
using namespace PDEs;
using namespace Physics;

namespace PDEs
{
  namespace Physics
  {
    CrossSections::CrossSections() :
        MaterialProperty(MaterialPropertyType::CROSS_SECTIONS)
    {}


    void CrossSections::reset()
    {
      n_groups = n_precursors = n_moments = 0;
      is_fissile = false; density = 1.0;

      e_bounds.clear();

      sigma_t.clear();
      sigma_a = sigma_s = sigma_f = sigma_r = sigma_t;
      chi = chi_prompt = sigma_t;
      nu = nu_prompt = nu_delayed = beta = sigma_t;
      nu_sigma_f = nu_prompt_sigma_f = nu_delayed_sigma_f = sigma_t;
      inv_velocity = diffusion_coeff = buckling = sigma_t;

      precursor_lambda.clear();
      precursor_yield.clear();
      chi_delayed.clear();

      transfer_matrices.clear();
    }


    void CrossSections::reinit()
    {
      assert(n_groups > 0);
      e_bounds.assign(n_groups + 1, 0.0);
      sigma_t.assign(n_groups, 0.0);
      sigma_a = sigma_s = sigma_f = sigma_r = sigma_t;
      chi = chi_prompt = sigma_t;
      nu = nu_prompt = nu_delayed = beta = sigma_t;
      nu_sigma_f = nu_prompt_sigma_f = nu_delayed_sigma_f = sigma_t;
      inv_velocity = diffusion_coeff = buckling = sigma_t;
    }


    void CrossSections::make_pure_scatterer()
    {
      assert(n_groups > 0);

      // set total cross-section to scattering cross-section
      sigma_t = sigma_s;

      // zero out capture reactions
      sigma_a.assign(sigma_a.size(), 0.0);
      sigma_r = sigma_f = sigma_a;
      chi = chi_prompt = sigma_a;
      nu = nu_prompt = nu_delayed = beta = sigma_a;
      nu_sigma_f = nu_prompt_sigma_f = nu_delayed_sigma_f = sigma_a;

      // recompute diffusion coefficient with standard definition
      for (unsigned int g = 0; g < n_groups; ++g)
        diffusion_coeff[g] = 1.0 / (3.0 * sigma_t[g]);

      is_fissile = false;
    }


    void CrossSections::compute_scattering_from_transfers()
    {
      sigma_s.assign(n_groups, 0.0);
      if (transfer_matrices.empty())
        return;

      // A group's scattering cross-section is defined as the sum of the
      // transfer cross-sections from a fixed group to all  other groups. In the
      // transfer matrices, the rows contain the destination groups and columns
      // the origin group. Due to this, computing sigma_s necessitates a column-
      // wise sum.
      for (unsigned int g = 0; g < n_groups; ++g)
        for (unsigned int gp = 0; gp < n_groups; ++gp)
        {
          assert(transfer_matrices[0][g][gp] >= 0.0);
          sigma_s[gp] += transfer_matrices[0][g][gp];
        }
    }


    void CrossSections::reconcile_cross_sections()
    {
      // determine whether sigma_a was specified
      double sum = accumulate(sigma_a.begin(), sigma_a.end(), 0.0);
      bool has_sigma_a = (sum > 1.0e-12);

      // Compute absorption xs from transfer matrix, if not specified. Otherwise,
      // recompute the total cross-section from the specified scattering and
      // absorption.
      if (not has_sigma_a)
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

      // compute the removal cross sections
      for (unsigned int g = 0; g < n_groups; ++g)
        sigma_r[g] = sigma_t[g] - transfer_matrices[0][g][g];
    }


    void CrossSections::reconcile_fission_properties()
    {
      //------------------------------------------------------------
      /// Lambda for finding non-zero cross-sections.
      auto is_nonzero_vector =
          [](std::vector<double> x)
          {
            auto is_nonzero = [](double y) { return y != 0.0; };
            return std::any_of(x.begin(), x.end(), is_nonzero);
          };

      // check which quantities were specified
      const bool has_sigf = is_nonzero_vector(sigma_f);
      const bool has_nusigf = is_nonzero_vector(nu_sigma_f);
      const bool has_nupsigf = is_nonzero_vector(nu_prompt_sigma_f);
      const bool has_nudsigf = is_nonzero_vector(nu_delayed_sigma_f);

      const bool has_nu = is_nonzero_vector(nu);
      const bool has_nup = is_nonzero_vector(nu_prompt);
      const bool has_nud = is_nonzero_vector(nu_delayed);
      const bool has_beta = is_nonzero_vector(beta);

      const bool has_chi = is_nonzero_vector(chi);
      const bool has_chip = is_nonzero_vector(chi_prompt);

      is_fissile = has_sigf || has_nusigf || has_nupsigf;

      // clear precursors if not fissile
      if (not is_fissile && n_precursors > 0)
      {
        n_precursors = 0;
        precursor_lambda.clear();
        precursor_yield.clear();
        chi_delayed.clear();

        std::cout << "!!! WARNING !!! " << "CrossSections::" << __FUNCTION__
                  << ": Precursors found in a non-fissile material.\n"
                  << "Clearing the precursor properties for consistency."
                  << std::endl;
      }

      // check fission properties
      if (is_fissile)
      {
        // Check the specified properties from most general to least.
        // Preference is given when the explicit properties are provided.
        // For example, if prompt and delayed quantities are explicitly
        // specified, the values are kept and others derived from them.
        // If explicit prompt and delayed  quantities are not specified,
        // alternate specification forms are checked. At present, the only
        // alternate way to specify prompt and delayed quantities is to use the
        // ''delayed fraction'' method where nu_sigma_f, nu, and beta are
        // specified group-wise. From these, all prompt and delayed quantities
        // can be derived. If neither of these conditions are met, it is
        // assumed that only total fission is desired and the prompt and
        // delayed properties are set accordingly.

        // begin with prompt/delayed quantities
        if (n_precursors > 0)
        {
          // prompt + delayed given
          if (has_sigf && has_nup && has_nud)
          {
            assert(all_of(sigma_f.begin(), sigma_f.end(),
                          [](double x) { return x >= 0.0; }));

            assert(all_of(nu_prompt.begin(), nu_prompt.end(),
                          [](double x) { return x > 0.0; }));

            assert(all_of(nu_delayed.begin(), nu_delayed.end(),
                          [](double x) { return x > 0.0; }));

            for (unsigned int g = 0; g < n_groups; ++g)
            {
              nu[g] = nu_prompt[g] + nu_delayed[g];
              beta[g] = nu_delayed[g] / nu[g];
            }
          }

          // delayed fraction given (1)
          else if (has_sigf && has_nu && has_beta)
          {
            assert(all_of(sigma_f.begin(), sigma_f.end(),
                          [](double x) { return x >= 0.0; }));

            assert(all_of(nu.begin(), nu.end(),
                          [](double x) { return x > 1.0; }));

            assert(all_of(beta.begin(), beta.end(),
                          [](double x) { return x > 0.0; }));

            for (unsigned int g = 0; g < n_groups; ++g)
            {
              nu_prompt[g] = (1.0 - beta[g]) * nu[g];
              nu_delayed[g] = beta[g] * nu[g];
            }
          }

          // delayed fraction given (2)
          else if (has_nusigf && has_nu && has_beta)
          {
            assert(all_of(nu_sigma_f.begin(), nu_sigma_f.end(),
                          [](double x) { return x >= 0.0; }));

            assert(all_of(nu.begin(), nu.end(),
                          [](double x) { return x >= 1.0; }));

            assert(all_of(beta.begin(), beta.end(),
                          [](double x) { return x > 0.0; }));

            for (unsigned int g = 0; g < n_groups; ++g)
            {
              sigma_f[g] = nu_sigma_f[g] / nu[g];
              nu_prompt[g] = (1.0 - beta[g]) * nu[g];
              nu_delayed[g] = beta[g] * nu[g];
            }
          }

          // no nu given
          if (has_nupsigf && has_nudsigf)
          {
            assert(all_of(nu_prompt_sigma_f.begin(), nu_prompt_sigma_f.end(),
                          [](double x) { return x >= 0.0; }));

            assert(all_of(nu_delayed_sigma_f.begin(), nu_delayed_sigma_f.end(),
                          [](double x) { return x >= 0.0; }));

            for (unsigned int g = 0; g < n_groups; ++g)
            {
              nu[g] = 1.0;
              sigma_f[g] = nu_prompt_sigma_f[g] + nu_delayed_sigma_f[g];
              nu_sigma_f[g] = sigma_f[g];

              beta[g] = nu_delayed_sigma_f[g] / nu_sigma_f[g];
              nu_prompt[g] = (1.0 - beta[g]);
              nu_delayed[g] = beta[g];
            }
          }

          // check the spectra
          assert(has_chip);

          assert(all_of(chi_prompt.begin(), chi_prompt.end(),
                        [](double x) { return x >= 0.0; }));

          assert(all_of(chi_delayed.begin(), chi_delayed.end(),
                        [](vector<double> x) {
                          return all_of(x.begin(), x.end(),
                                        [](double y) { return y >= 0.0; });
                        }));

          // check precursor properties
          assert(all_of(precursor_yield.begin(), precursor_yield.end(),
                        [](double x) { return x >= 0.0; }));

          assert(all_of(precursor_lambda.begin(), precursor_lambda.end(),
                        [](double x) { return x > 0.0; }));

          // normalize prompt spectrum
          double prompt_sum = accumulate(
              chi_prompt.begin(), chi_prompt.end(), 0.0
          );

          if (abs(prompt_sum - 1.0) > 1.0e-12)
          {
            stringstream warn;
            warn << "**!!!** Warning **!!!** CrossSections::"
                 << __FUNCTION__ << ": "
                 << "Normalizing prompt fission spectrum to unity.";
            std::cout << warn.str() << endl;

            for (unsigned int g = 0; g < n_groups; ++g)
              chi_prompt[g] /= prompt_sum;
          }

          // normalize delayed spectra
          for (unsigned int j = 0; j < n_precursors; ++j)
          {
            double delayed_sum_j = 0.0;
            for (unsigned int g = 0; g < n_groups; ++g)
              delayed_sum_j += chi_delayed[g][j];

            if (abs(delayed_sum_j - 1.0) > 1.0e-12)
            {
              stringstream warn;
              warn << "**!!!** Warning **!!!** CrossSections::"
                   << __FUNCTION__ << ": "
                   << "Normalizing delayed emission spectrum to unity "
                   << "for precursor species " << j << ".";
              std::cout << warn.str() << endl;
            }

            for (unsigned int g = 0; g < n_groups; ++g)
              chi_delayed[g][j] /= delayed_sum_j;
          }

          // normalize precursor yeilds
          double yield_sum = accumulate(
              precursor_yield.begin(), precursor_yield.end(), 0.0
          );

          if (abs(yield_sum - 1.0) > 1.0e-12)
          {
            stringstream warn;
            warn << "**!!!** Warning **!!!** CrossSections::"
                 << __FUNCTION__ << ": "
                 << "Normalizing precursor yields to unity.";
            std::cout << warn.str() << endl;

            for (unsigned int j = 0; j < n_precursors; ++j)
              precursor_yield[j] /= yield_sum;
          }

          // compute total quantities
          for (unsigned int g = 0; g < n_groups; ++g)
          {
            nu[g] = nu_prompt[g] + nu_delayed[g];

            // compute the beta-weighted total fission spectra
            chi[g] = (1.0 - beta[g]) * chi_prompt[g];
            for (unsigned int j = 0; j < n_precursors; ++j)
              chi[g] += beta[g] * precursor_yield[j] * chi_delayed[g][j];
          }

          // check results of above
          assert(all_of(nu.begin(), nu.end(),
                        [](double x) { return x >= 1.0; }));

          assert(all_of(chi.begin(), chi.end(),
                        [](double x) { return x >= 0.0; }));

          double chi_sum = accumulate(chi.begin(), chi.end(), 0.0);
          for (unsigned int g = 0; g < n_groups; ++g)
            chi[g] /= chi_sum;
        }//if precursors

          // total fission only
        else
        {
          assert(has_chi && has_nu && (has_nusigf || has_sigf));

          assert(all_of(nu.begin(), nu.end(),
                        [](double x) { return x >= 1.0; }));

          assert(all_of(chi.begin(), chi.end(),
                        [](double x) { return x >= 0.0; }));

          if (has_nusigf && has_nu && not has_sigf)
            for (unsigned int g = 0; g < n_groups; ++g)
              sigma_f[g] = nu_sigma_f[g] / nu[g];

          assert(all_of(sigma_f.begin(), sigma_f.end(),
                        [](double x) { return x >= 0.0; }));

          // normalize the fission spectrum
          double sum = accumulate(chi.begin(), chi.end(), 0.0);
          if (abs(sum - 1.0) > 1.0e-12)
          {
            stringstream warn;
            warn << "**!!!** Warning **!!!** CrossSections::"
                 << __FUNCTION__ << ": "
                 << "Normalizing total fission spectrum to unity.";
            std::cout << warn.str() << endl;

            for (unsigned int g = 0; g < n_groups; ++g)
              chi[g] /= sum;
          }
        }

        // compute nu_sigma_f terms
        for (unsigned int g = 0; g < n_groups; ++g)
        {
          nu_sigma_f[g] = nu[g] * sigma_f[g];
          nu_prompt_sigma_f[g] = nu_prompt[g] * sigma_f[g];
          nu_delayed_sigma_f[g] = nu_delayed[g] * sigma_f[g];
        }
      }//if fissile
    }


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

      for (unsigned int ell = 0; ell < n_moments; ++ell)
        for (unsigned int g = 0; g < n_groups; ++g)
          for (unsigned int gp = 0; gp < n_groups; ++gp)
            transfer_matrices[ell][g][gp] *= density;

      // compute diffusion coefficient if unspecified
      double sum = accumulate(diffusion_coeff.begin(),
                              diffusion_coeff.end(), 0.0);
      if (sum < 1.0e-12)
      {
        for (unsigned int g = 0; g < n_groups; ++g)
          diffusion_coeff[g] = 1.0 / (3.0 * sigma_t[g]);
      }
    }

  }
}
