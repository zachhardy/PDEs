#include "transient_solver.h"


using namespace NeutronDiffusion;


void
TransientSolver::set_transient_source(Groupset& groupset,
                                      SourceFlags source_flags)
{
  if (source_flags == NO_SOURCE_FLAGS)
    return;

  const bool apply_mat_src = (source_flags & APPLY_MATERIAL_SOURCE);
  const bool apply_wgs_scatter_src = (source_flags & APPLY_WGS_SCATTER_SOURCE);
  const bool apply_ags_scatter_src = (source_flags & APPLY_AGS_SCATTER_SOURCE);
  const bool apply_wgs_fission_src = (source_flags & APPLY_WGS_FISSION_SOURCE);
  const bool apply_ags_fission_src = (source_flags & APPLY_AGS_FISSION_SOURCE);

  // Get groupset range
  const auto n_gsg = groupset.groups.size();
  const auto gs_i = groupset.groups.front();
  const auto gs_f = groupset.groups.back();

  // Get timestep size
  const double eff_dt = effective_time_step();

  // Get groupset vector
  auto& b = groupset.b;

  //======================================== Loop over cells
  for (const auto& cell : mesh->cells)
  {
    const double volume = cell.volume;
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];

    const size_t uk_map = cell.id * n_groups;
    const size_t gs_uk_map = cell.id * n_gsg;
    const size_t prec_uk_map = cell.id * max_precursors;

    const auto src_id = matid_to_src_map[cell.material_id];
    const double* src = nullptr;
    if (src_id >= 0 && apply_mat_src)
      src = material_src[src_id]->values.data();

    //============================== Loop over groups
    for (size_t gr = 0; gr < n_gsg; ++gr)
    {
      double rhs = 0.0;
      const size_t g = groupset.groups[gr];

      const double* inv_vel = xs->inv_velocity.data();

      //==================== Inhomogeneous source term
      rhs += (src)? src[g] : 0.0;

      //==================== Old scalar flux
      rhs += inv_vel[g]/eff_dt * phi_old[uk_map + g];

      //==================== Old precursors
      if (xs->is_fissile && use_precursors)
      {
        const double* chi_d = xs->chi_delayed[g].data();
        const double* lambda = xs->precursor_lambda.data();
        for (size_t j = 0; j < xs->n_precursors; ++j)
        {
          double coeff = chi_d[j] * lambda[j];
          if (not lag_precursors)
            coeff /= 1.0 + eff_dt * lambda[j];

          rhs += coeff * precursor_old[prec_uk_map + j];
        }
      }//if precursors

      //==================== Within-groupset scattering
      if (apply_wgs_scatter_src)
      {
        const double* sig_s = xs->transfer_matrices[0][g].data();

        for (const auto& gp : groupset.groups)
          rhs += sig_s[gp] * phi[uk_map + gp];
      }//if within-groupset scattering

      //==================== Across-groupset scattering
      if (apply_ags_scatter_src)
      {
        const double* sig_s = xs->transfer_matrices[0][g].data();

        for (size_t gpr = 0; gpr < n_groups; ++gpr)
        {
          const size_t gp = groups[gpr];
          if (gp < gs_i || gp > gs_f)
            rhs += sig_s[gp] * phi[uk_map + gp];
        }//for gprime
      }//if across-groupset scattering

      //==================== Within-groupset fission
      if (xs->is_fissile && apply_wgs_fission_src)
      {
        //==================== Total fission
        if (not use_precursors)
        {
          const double chi = xs->chi[g];
          const double* nu_sigf = xs->nu_sigma_f.data();

          for (const auto& gp : groupset.groups)
            rhs += chi * nu_sigf[gp] * phi[uk_map + gp];
        }//if total fission

        //==================== Prompt + delayed fission
        else
        {
          const double chi_p = xs->chi_prompt[g];
          const double* chi_d = xs->chi_delayed[g].data();
          const double* nup_sigf = xs->nu_prompt_sigma_f.data();
          const double* nud_sigf = xs->nu_delayed_sigma_f.data();
          const double* lambda = xs->precursor_lambda.data();
          const double* gamma = xs->precursor_yield.data();

          //===== Prompt fission
          for (const auto& gp : groupset.groups)
            rhs += chi_p * nup_sigf[gp] * phi[uk_map + gp];

          //===== Delayed fission
          if (not lag_precursors)
          {
            double coeff = 0.0;
            for (size_t j = 0; j < xs->n_precursors; ++j)
              coeff += chi_d[j] * lambda[j] / (1.0 + eff_dt*lambda[j]) *
                       gamma[j] * eff_dt;

            for (const auto& gp : groupset.groups)
              rhs += coeff * nud_sigf[gp] * phi[uk_map + gp];
          }//if not lagging precursors
        }//if prompt+delayed fission
      }//if within-groupset fission

      if (xs->is_fissile && apply_ags_fission_src)
      {
        //==================== Total fission
        if (not use_precursors)
        {
          const double chi = xs->chi[g];
          const double* nu_sigf = xs->nu_sigma_f.data();

          for (size_t gpr = 0; gpr < n_groups; ++gpr)
          {
            const size_t gp = groups[gpr];
            if (gp < gs_i || gp > gs_f)
              rhs += chi * nu_sigf[gp] * phi[uk_map + gp];
          }//for gprime
        }//if total fission

          //==================== Prompt + delayed fission
        else
        {
          const double chi_p = xs->chi_prompt[g];
          const double* chi_d = xs->chi_delayed[g].data();
          const double* nup_sigf = xs->nu_prompt_sigma_f.data();
          const double* nud_sigf = xs->nu_delayed_sigma_f.data();
          const double* lambda = xs->precursor_lambda.data();
          const double* gamma = xs->precursor_yield.data();

          //===== Prompt fission
          for (size_t gpr = 0; gpr < n_groups; ++gpr)
          {
            const size_t gp = groups[gpr];
            if (gp < gs_i || gp > gs_f)
              rhs += chi_p * nup_sigf[gp] * phi[uk_map + gp];
          }//for gprime

          //===== Delayed fission
          if (not lag_precursors)
          {
            double coeff = 0.0;
            for (size_t j = 0; j < xs->n_precursors; ++j)
              coeff += chi_d[j] * lambda[j] / (1.0 + eff_dt*lambda[j]) *
                       gamma[j] * eff_dt;

            for (size_t gpr = 0; gpr < n_groups; ++gpr)
            {
              const size_t gp = groups[gpr];
              if (gp < gs_i || gp > gs_f)
                rhs += coeff * nud_sigf[gp] * phi[uk_map + gp];
            }//for gprime
          }//if not lagging precursors
        }//if prompt+delayed fission
      }//if across-groupset fission

      b[gs_uk_map + gr] += rhs * volume;

    }//for group

    //======================================== Loop over faces
    for (const auto& face : cell.faces)
    {
      // Skip interior faces
      if (face.has_neighbor)
        continue;

      const auto bndry_id = face.neighbor_id;
      const auto& bndry_type = boundary_info[bndry_id].first;

      //==================== Dirichlet term
      if (bndry_type == BoundaryType::DIRICHLET)
      {
        const double* D = xs->diffusion_coeff.data();
        const double d_pf = cell.centroid.distance(face.centroid);

        for (size_t gr = 0; gr < n_gsg; ++gr)
        {
          const size_t g = groupset.groups[gr];

          const auto& bndry = boundaries[bndry_id][gr];
          const auto bc = std::static_pointer_cast<DirichletBoundary>(bndry);

          b[gs_uk_map + gr] += D[g]/d_pf * bc->value * face.area;
        }
      }//if Dirichlet

      //==================== Neumann term
      if (bndry_type == BoundaryType::NEUMANN)
      {
        for (size_t gr = 0; gr < n_gsg; ++gr)
        {
          const size_t g = groupset.groups[gr];

          const auto& bndry = boundaries[bndry_id][gr];
          const auto bc = std::static_pointer_cast<NeumannBoundary>(bndry);

          b[gs_uk_map + gr] += bc->value*face.area;
        }
      }//if Neumann


      //==================== Robin term
      else if (bndry_type == BoundaryType::MARSHAK ||
               bndry_type == BoundaryType::ROBIN)
      {
        const double* D = &xs->diffusion_coeff[0];
        const double d_pf = cell.centroid.distance(face.centroid);

        for (size_t gr = 0; gr < n_gsg; ++gr)
        {
          const size_t g = groupset.groups[gr];

          const auto& bndry = boundaries[bndry_id][gr];
          const auto bc = std::static_pointer_cast<RobinBoundary>(bndry);

          const double coeff = D[g]/(bc->b*D[g] + bc->a*d_pf);
          b[gs_uk_map + gr] += coeff * bc->f * face.area;
        }
      }//if Robin
    }//for face
  }//for cell
}
