#include "../steadystate_solver.h"


using namespace NeutronDiffusion;


void
SteadyStateSolver::
fv_set_source(Groupset& groupset, SourceFlags source_flags)
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

  // Get groupset vector
  auto& b = groupset.b;

  //================================================== Loop over cells
  for (const auto& cell : mesh->cells)
  {
    const double volume = cell.volume;
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];

    const size_t i = cell.id*n_gsg;
    const size_t uk_map = cell.id*n_groups;

    const double* x = &phi[uk_map];

    //======================================== Inhomogeneous source term
    const int src_id = matid_to_src_map[cell.material_id];
    const double* src = nullptr;
    if (src_id >= 0 and apply_mat_src)
      src = &material_src[src_id]->values[0];

    //======================================== Group coupling terms
    for (size_t gr = 0; gr < n_gsg; ++gr)
    {
      const size_t ig = i + gr;
      const size_t g = groupset.groups[gr];

      double rhs = 0.0;

      //==================== Inhomogeneous source term
      if (src) rhs += src[g]*volume;

      //==================== Scattering terms
      if (apply_wgs_scatter_src || apply_ags_scatter_src)
      {
        const double* sig_s = &xs->transfer_matrices[0][gr][0];

        if (apply_wgs_scatter_src)
          for (const auto& gp : groupset.groups)
            rhs += sig_s[gp]*x[gp]*volume;

        if (apply_ags_scatter_src)
          for (size_t gpr = 0; gpr < n_groups; ++gpr)
          {
            const size_t gp = groups[gpr];
            if (gp < gs_i || gp > gs_f)
              rhs += sig_s[gp]*x[gp]*volume;
          }
      }

      //==================== Fission terms
      if (xs->is_fissile && apply_wgs_fission_src || apply_ags_fission_src)
      {
        //==================== Total fission
        if (not use_precursors)
        {
          const double chi = xs->chi[g];
          const double* nu_sigf = &xs->nu_sigma_f[0];

          if (apply_wgs_fission_src)
            for (const auto& gp : groupset.groups)
              rhs += chi*nu_sigf[gp]*x[gp]*volume;

          if (apply_ags_fission_src)
            for (size_t gpr = 0; gpr < n_groups; ++gpr)
            {
              const size_t gp = groups[gpr];
              if (gp < gs_i || gp > gs_f)
                rhs += chi*nu_sigf[gp]*x[gp]*volume;
            }
        }

          //==================== Prompt + delayed fission
        else
        {
          const double chi_p = xs->chi_prompt[g];
          const double* chi_d = &xs->chi_delayed[g][0];
          const double* nup_sigf = &xs->nu_prompt_sigma_f[0];
          const double* nud_sigf = &xs->nu_delayed_sigma_f[0];
          const double* gamma = &xs->precursor_yield[0];

          if (apply_wgs_fission_src)
          {
            for (const auto& gp : groupset.groups)
            {
              double value = chi_p*nup_sigf[gp];
              for (size_t j = 0; j < xs->n_precursors; ++j)
                value += chi_d[j]*gamma[j]*nud_sigf[gp];
              rhs += value*x[gp]*volume;
            }
          }

          if (apply_ags_fission_src)
          {
            for (size_t gpr = 0; gpr < n_groups; ++gpr)
            {
              const size_t gp = groups[gpr];
              if (gp < gs_i || gp > gs_f)
              {
                double value = chi_p*nup_sigf[gp];
                for (size_t j = 0; j < xs->n_precursors; ++j)
                  value += chi_d[j]*gamma[j]*nud_sigf[gp];
                rhs += value*x[gp]*volume;
              }
            }//for gpr
          }
        }//prompt + delayed fission
      }//if fissile

      b[ig] += rhs;

    }//for gr

    //================================================== Loop over faces
    for (const auto& face : cell.faces)
    {
      //============================== Boundary faces
      if (not face.neighbor_id)
      {
        const auto bndry_id = face.neighbor_id;
        const auto& bndry_type = boundary_info[bndry_id].first;

        //==================== Dirichlet boundary term
        if (bndry_type == BoundaryType::DIRICHLET)
        {
          const double* D = &xs->diffusion_coeff[0];
          const double d_pf = cell.centroid.distance(face.centroid);

          for (size_t gr = 0; gr < n_gsg; ++gr)
          {
            const size_t ig = i + gr;
            const size_t g = groupset.groups[gr];

            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<DirichletBoundary>(bndry);

            b[ig] += D[g]/d_pf*bc->value*face.area;
          }
        }

          //==================== Neumann boundary term
        else if (bndry_type == BoundaryType::NEUMANN)
        {
          for (size_t gr = 0; gr < n_gsg; ++gr)
          {
            const size_t ig = i + gr;
            const size_t g = groupset.groups[gr];

            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<NeumannBoundary>(bndry);

            b[ig] += bc->value*face.area;
          }
        }

          //==================== Robin boundary term
        else if (bndry_type == BoundaryType::MARSHAK or
                 bndry_type == BoundaryType::ROBIN)
        {
          const double* D = &xs->diffusion_coeff[0];
          const double d_pf = cell.centroid.distance(face.centroid);

          for (size_t gr = 0; gr < n_gsg; ++gr)
          {
            const size_t ig = i + gr;
            const size_t g = groupset.groups[gr];

            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<RobinBoundary>(bndry);

            b[ig] += D[g]/(bc->b*D[g] + bc->a*d_pf)*bc->f*face.area;
          }
        }
      }//if boundary face
    }//for faces
  }//for cells
}
