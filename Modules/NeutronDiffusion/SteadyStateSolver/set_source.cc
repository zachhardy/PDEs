#include "steadystate_solver.h"


using namespace NeutronDiffusion;


void
SteadyStateSolver::
set_source(SourceFlags source_flags)
{
  if (source_flags == NO_SOURCE_FLAGS)
    return;

  const bool apply_mat_src = (source_flags & APPLY_MATERIAL_SOURCE);
  const bool apply_scatter_src = (source_flags & APPLY_SCATTER_SOURCE);
  const bool apply_fission_src = (source_flags & APPLY_FISSION_SOURCE);
  const bool apply_bndry_src = (source_flags & APPLY_BOUNDARY_SOURCE);

  // Loop over cells
  for (const auto& cell : mesh->cells)
  {
    const auto volume = cell.volume;
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];

    size_t uk_map = n_groups * cell.id;

    const auto src_id = matid_to_src_map[cell.material_id];
    const double* src = nullptr;
    if (src_id >= 0 && apply_mat_src)
      src = material_src[src_id]->values.data();

    // Loop over groups
    for (unsigned int gr = 0; gr < n_groups; ++gr)
    {
      double rhs = 0.0;
      const auto g = groups[gr];

      //========================================
      // Inhomogeneous source term
      //========================================

      rhs += (src)? src[g] : 0.0;

      //========================================
      // Scattering source term
      //========================================

      if (apply_scatter_src)
      {
        const auto* sig_s = xs->transfer_matrices[0][g].data();
        for (unsigned int gpr = 0; gpr < n_groups; ++gpr)
          rhs += sig_s[groups[gpr]] * phi[uk_map + gpr];
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
          for (unsigned int gpr = 0; gpr < n_groups; ++gpr)
            rhs += chi * nu_sigf[groups[gpr]] * phi[uk_map + gpr];
        }//if total fission

        // Prompt + delayed fission
        else
        {
          const auto chi_p = xs->chi_prompt[g];
          const auto* chi_d = xs->chi_delayed[g].data();
          const auto* nup_sigf = xs->nu_prompt_sigma_f.data();
          const auto* nud_sigf = xs->nu_delayed_sigma_f.data();
          const auto* gamma = xs->precursor_yield.data();

          // Prompt fission
          for (unsigned int gpr = 0; gpr < n_groups; ++gpr)
            rhs += chi_p * nup_sigf[groups[gpr]] * phi[uk_map + gpr];

          // Delayed fission
          double coeff = 0.0;
          for (unsigned int j = 0; j < xs->n_precursors; ++j)
            coeff += chi_d[j] * gamma[j];

          for (unsigned int gpr = 0; gpr < n_groups; ++gpr)
            rhs += coeff * nud_sigf[groups[gpr]] * phi[uk_map + gpr];
        }//if prompt+delayed fission
      }//if fission

      b[uk_map + gr] += rhs * volume;

    }//for group

    //========================================
    // Boundary source terms
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

          for (unsigned int gr = 0; gr < n_groups; ++gr)
          {
            const auto &bndry = boundaries[bndry_id][gr];
            const auto bc = std::static_pointer_cast<DirichletBoundary>(bndry);

            b[uk_map + gr] += D[groups[gr]] / d_pf * bc->value * face.area;
          }
        }//if Dirichlet

        //========================================
        // Neumann boundary term
        //========================================

        if (bndry_type == BoundaryType::NEUMANN)
        {
          for (unsigned int gr = 0; gr < n_groups; ++gr)
          {
            const auto &bndry = boundaries[bndry_id][gr];
            const auto bc = std::static_pointer_cast<NeumannBoundary>(bndry);

            b[uk_map + gr] += bc->value * face.area;
          }
        }//if Neumann

          //========================================
          // Robin boundary terms
          //========================================

        else if (bndry_type == BoundaryType::MARSHAK ||
                 bndry_type == BoundaryType::ROBIN)
        {
          const auto *D = xs->diffusion_coeff.data();
          const auto d_pf = cell.centroid.distance(face.centroid);

          for (unsigned int gr = 0; gr < n_groups; ++gr)
          {
            const auto g = groups[gr];
            const auto &bndry = boundaries[bndry_id][gr];
            const auto bc = std::static_pointer_cast<RobinBoundary>(bndry);

            const double coeff = D[g] / (bc->b * D[g] + bc->a * d_pf);
            b[uk_map + gr] += coeff * bc->f * face.area;
          }
        }//if Robin
      }//for face
    }//if boundary source
  }//for cell
}
