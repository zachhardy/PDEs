#include "transient_solver.h"


using namespace NeutronDiffusion;


void
TransientSolver::
set_transient_source(SourceFlags source_flags)
{
  if (source_flags == NO_SOURCE_FLAGS)
    return;

  const bool apply_mat_src = (source_flags & APPLY_MATERIAL_SOURCE);
  const bool apply_scatter_src = (source_flags & APPLY_SCATTER_SOURCE);
  const bool apply_fission_src = (source_flags & APPLY_FISSION_SOURCE);
  const bool apply_bndry_src = (source_flags & APPLY_BOUNDARY_SOURCE);

  // Get timestep size
  const double eff_dt = effective_time_step();

  // Loop over cells
  for (const auto& cell : mesh->cells)
  {
    const auto volume = cell.volume;
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];

    const auto uk_map_g = n_groups * cell.id;
    const auto uk_map_j = max_precursors * cell.id;

    const auto src_id = matid_to_src_map[cell.material_id];
    const double* src = nullptr;
    if (src_id < material_src.size() && apply_mat_src)
      src = material_src[src_id]->values.data();

    // Loop over groups
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      double rhs = 0.0;
      const auto* inv_vel = xs->inv_velocity.data();

      //========================================
      // Inhomogeneous source term
      //========================================

      rhs += (src)? src[g] : 0.0;

      //========================================
      // Old scalar flux
      //========================================

      rhs += inv_vel[g]/eff_dt * phi_old[uk_map_g + g];

      //========================================
      // Old precursors
      //========================================

      if (xs->is_fissile && use_precursors)
      {
        const auto* chi_d = xs->chi_delayed[g].data();
        const auto* lambda = xs->precursor_lambda.data();
        for (unsigned int j = 0; j < xs->n_precursors; ++j)
        {
          double coeff = chi_d[j] * lambda[j];
          if (not lag_precursors)
            coeff /= 1.0 + eff_dt * lambda[j];

          rhs += coeff*precursors_old[uk_map_j + j];
        }
      }//if precursors

      //========================================
      // Scattering source term
      //========================================
      if (apply_scatter_src)
      {
        const auto* sig_s = xs->transfer_matrices[0][g].data();
        for (unsigned int gp = 0; gp < n_groups; ++gp)
          rhs += sig_s[gp] * phi[uk_map_g + gp];
      }//if scattering

      //========================================
      // Fission source term
      //========================================

      if (xs->is_fissile && apply_fission_src)
      {
        // Total fission
        if (not use_precursors)
        {
          const auto chi = xs->chi[g];
          const auto* nu_sigf = xs->nu_sigma_f.data();
          for (unsigned int gp = 0; gp < n_groups; ++gp)
            rhs += chi * nu_sigf[gp] * phi[uk_map_g + gp];
        }//if total fission

        // Prompt + delayed fission
        else
        {
          const auto chi_p = xs->chi_prompt[g];
          const auto* chi_d = xs->chi_delayed[g].data();
          const auto* nup_sigf = xs->nu_prompt_sigma_f.data();
          const auto* nud_sigf = xs->nu_delayed_sigma_f.data();
          const auto* lambda = xs->precursor_lambda.data();
          const auto* gamma = xs->precursor_yield.data();

          // Prompt
          for (unsigned int gp = 0; gp < n_groups; ++gp)
            rhs += chi_p * nup_sigf[gp] * phi[uk_map_g + gp];

          // Delayed
          if (not lag_precursors)
          {
            double coeff = 0.0;
            for (unsigned int j = 0; j < xs->n_precursors; ++j)
              coeff += chi_d[j] * lambda[j] / (1.0 + eff_dt*lambda[j]) *
                       gamma[j] * eff_dt;

            for (unsigned int gp = 0; gp < n_groups; ++gp)
              rhs += coeff * nud_sigf[gp] * phi[uk_map_g + gp];
          }//if not lagging precursors
        }//if prompt+delayed fission
      }//if fission

      b[uk_map_g + g] += rhs * volume;

    }//for group

    //========================================
    // Boundary source term
    //========================================
    if (apply_bndry_src)
    {
      // Loop over faces
      for (const auto &face: cell.faces)
      {
        // Skip interior faces
        if (face.has_neighbor)
          continue;

        const auto bndry_id = face.neighbor_id;
        const auto &bndry_type = boundary_info[bndry_id].first;

        //========================================
        // Dirichlet boundary term
        //========================================

        if (bndry_type == BoundaryType::DIRICHLET)
        {
          const auto *D = xs->diffusion_coeff.data();
          const auto d_pf = cell.centroid.distance(face.centroid);

          for (unsigned int g = 0; g < n_groups; ++g)
          {
            const auto &bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<DirichletBoundary>(bndry);
            b[uk_map_g + g] += D[g] / d_pf * bc->value * face.area;
          }
        }//if Dirichlet

        //========================================
        // Neumann boundary term
        //========================================

        if (bndry_type == BoundaryType::NEUMANN)
        {
          for (unsigned int g = 0; g < n_groups; ++g)
          {
            const auto &bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<NeumannBoundary>(bndry);
            b[uk_map_g + g] += bc->value * face.area;
          }
        }//if Neumann

          //========================================
          // Robin boundary term
          //========================================

        else if (bndry_type == BoundaryType::MARSHAK ||
                 bndry_type == BoundaryType::ROBIN)
        {
          const auto *D = xs->diffusion_coeff.data();
          const auto d_pf = cell.centroid.distance(face.centroid);
          for (unsigned int g = 0; g < n_groups; ++g)
          {
            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<RobinBoundary>(bndry);

            const double coeff = D[g] / (bc->b * D[g] + bc->a * d_pf);
            b[uk_map_g + g] += coeff * bc->f * face.area;
          }
        }//if Robin
      }//for face
    }
  }//for cell
}
