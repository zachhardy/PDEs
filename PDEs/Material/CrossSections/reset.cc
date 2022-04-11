#include "cross_sections.h"

//##################################################
void CrossSections::reset()
{
  n_groups = 0;
  n_precursors = 0;
  scattering_order = 0;

  is_fissile = false;
  has_precursors = false;
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
