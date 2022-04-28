#include "steadystate_solver_fv.h"

/**
* \brief Assemble the matrix for the specified \p groupset.
 *
 * If solving the full system, this routine assembles the full multigroup
 * operator for the groupset, including off-diagonal scattering and fission
 * coupling terms. Otherwise, this routine assembles the within-group system
 * for the groups within the groupset.
*/
void neutron_diffusion::SteadyStateSolver_FV::
assemble_matrix(Groupset& groupset)
{
  math::Matrix<double>& A = groupset.matrix;
  A *= 0.0;

  // Get groupset range
  const auto n_gsg = groupset.groups.size();
  const auto gs_i = groupset.groups.front();
  const auto gs_f = groupset.groups.back();

  //================================================== Loop over cells
  for (const auto& cell : mesh->cells)
  {
    const double volume = cell->volume;
    const auto& xs = material_xs[matid_to_xs_map[cell->material_id]];
    const size_t i = cell->id * n_gsg;

    //==================== Total interaction terms
    for (size_t g = gs_i; g <= gs_f; ++g)
      A[i + g][i + g] += xs->sigma_t[g] * volume;

    //============================== Cross-group coupling
    if (solution_technique == SolutionTechnique::FULL_SYSTEM)
    {
      for (size_t g = gs_i; g <= gs_f; ++g)
      {
        //==================== Scattering term
        for (size_t gp = gs_i; gp <= gs_f; ++gp)
          A[i + g][i + gp] -= xs->transfer_matrices[0][g][gp] * volume;

        //==================== Fission term
        if (xs->is_fissile)
        {
          //========== Total fission
          if (not use_precursors)
          {
            const double chi = xs->chi[g];
            for (size_t gp = gs_i; gp <= gs_f; ++gp)
              A[i + g][i + gp] -= chi * xs->nu_sigma_f[gp] * volume;
          }

          //========== Prompt + delayed fission
          else
          {
            const double chi_p = xs->chi_prompt[g];
            for (size_t gp = gs_i; gp <= gs_f; ++gp)
              A[i + g][i + gp] -= chi_p * xs->nu_prompt_sigma_f[gp] * volume;

            for (size_t j = 0; j < xs->n_precursors; ++j)
            {
              const double chi_d = xs->chi_delayed[g][j];
              const double gamma = xs->precursor_yield[j];
              for (size_t gp = gs_i; gp <= gs_f; ++gp)
                A[i + g][i + gp] -= chi_d * gamma *
                                    xs->nu_delayed_sigma_f[gp] *
                                    volume;
            }//for precursors
          }
        }//if fissile
      }//for group
    }//if full system

    //================================================== Loop over faces
    for (const auto& face : cell->faces)
    {
      //============================== Interior faces
      if (face.has_neighbor)
      {
        // Get neighbor info
        const auto& nbr_cell = mesh->cells[face.neighbor_id];
        const size_t j = nbr_cell->id * n_gsg;

        const int nbr_xs_id = matid_to_xs_map[nbr_cell->material_id];
        const auto& nbr_xs = material_xs[nbr_xs_id];

        // Geometric quantities
        const double d_pf = cell->centroid.distance(face.centroid);
        const double d_pn = cell->centroid.distance(nbr_cell->centroid);
        const double w = d_pf / d_pn; // harmonic mean weighting factor

        //==================== Diffusion term
        for (size_t g = gs_i; g <= gs_f; ++g)
        {
          const double D = xs->diffusion_coeff[g];
          const double D_nbr = nbr_xs->diffusion_coeff[g];
          const double D_eff = 1.0 / (w/D + (1.0 - w)/D_nbr);

          A[i + g][i + g] += D_eff / d_pn * face.area;
          A[i + g][j + g] -= D_eff / d_pn * face.area;
        }
      }//if interior face

      //============================== Boundary faces
      else
      {
        const auto bndry_id = face.neighbor_id;
        const auto bndry_type = boundary_info[bndry_id].first;

        //==================== Dirichlet boundary term
        if (bndry_type == BoundaryType::ZERO_FLUX or
            bndry_type == BoundaryType::DIRICHLET)
        {
          const double d_pf = cell->centroid.distance(face.centroid);
          for (size_t g = gs_i; g <= gs_f; ++g)
          {
            const double D = xs->diffusion_coeff[g];
            A[i + g][i + g] += D / d_pf * face.area;
          }
        }

        //==================== Robin boundary term
        else if (bndry_type == BoundaryType::VACUUM or
                 bndry_type == BoundaryType::MARSHAK or
                 bndry_type == BoundaryType::ROBIN)
        {
          const double d_pf = cell->centroid.distance(face.centroid);
          for (size_t g = gs_i; g <= gs_f; ++g)
          {
            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<RobinBoundary>(bndry);
            const double D = xs->diffusion_coeff[g];
            A[i + g][i + g] += bc->a*D/(bc->b*D + bc->a*d_pf) * face.area;
          }
        }
      }//if boundary face
    }//for face
  }//for cell
}

//######################################################################

void neutron_diffusion::SteadyStateSolver_FV::
set_source(Groupset& groupset, math::Vector<double>& b,
           SourceFlags source_flags)
{
  const bool apply_mat_src         = (source_flags & APPLY_MATERIAL_SOURCE);
  const bool apply_wgs_scatter_src = (source_flags & APPLY_WGS_SCATTER_SOURCE);
  const bool apply_ags_scatter_src = (source_flags & APPLY_AGS_SCATTER_SOURCE);
  const bool apply_wgs_fission_src = (source_flags & APPLY_WGS_FISSION_SOURCE);
  const bool apply_ags_fission_src = (source_flags & APPLY_AGS_FISSION_SOURCE);

  // Get groupset range
  const auto n_gsg = groupset.groups.size();
  const auto gs_i = groupset.groups.front();
  const auto gs_f = groupset.groups.back();

  const auto g_i = groups.front();
  const auto g_f = groups.back();

  //================================================== Loop over cells
  for (const auto& cell : mesh->cells)
  {
    const double volume = cell->volume;
    const auto& xs = material_xs[matid_to_xs_map[cell->material_id]];

    const size_t i = cell->id * n_gsg;
    const size_t uk_map = cell->id * n_groups;

    //==================== Inhomogeneous source term
    const int src_id = matid_to_src_map[cell->material_id];
    if (src_id >= 0 and apply_mat_src)
    {
      const auto& src = material_src[src_id];
      for (size_t g = gs_i; g <= gs_f; ++g)
        b[i + g] += src->values[g] * volume;
    }

    for (size_t g = gs_i; g <= gs_f; ++g)
    {
      //============================== Within groupset scattering term
      if (apply_wgs_scatter_src)
        for (size_t gp = gs_i; g <= gs_f; ++g)
          b[i + g] += xs->transfer_matrices[0][g][gp] *
                      phi[uk_map + gp] * volume;

      //============================== Across groupset scattering term
      if (apply_ags_scatter_src)
        for (size_t gp = g_i; gp <= g_f; ++gp)
          if (gp < gs_i or gp > gs_f)
            b[i + g] += xs->transfer_matrices[0][g][gp] *
                        phi[uk_map + gp] * volume;

      //============================== Within groupset fission term
      if (xs->is_fissile and apply_wgs_fission_src)
      {
        if (not use_precursors)
        {
          const double chi = xs->chi[g];
          for (size_t gp = gs_i; gp <= gs_f; ++gp)
            b[i + g] += chi * xs->nu_sigma_f[gp] *
                        phi[uk_map + gp] * volume;
        }
        else
        {
          const double chi_p = xs->chi_prompt[g];
          for (size_t gp = gs_i; gp <= gs_f; ++gp)
            b[i + g] += chi_p * xs->nu_prompt_sigma_f[gp] *
                        phi[uk_map + gp] * volume;

          for (size_t j = 0; j < xs->n_precursors; ++j)
          {
            const double chi_d = xs->chi_delayed[g][j];
            const double gamma = xs->precursor_yield[j];
            for (size_t gp = gs_i; gp < gs_f; ++gp)
              b[i + g] += chi_d * gamma *
                          xs->nu_delayed_sigma_f[gp] *
                          phi[uk_map + gp] * volume;
          }
        }
      }

      if (xs->is_fissile and apply_wgs_fission_src)
      {
        if (not use_precursors)
        {
          const double chi = xs->chi[g];
          for (size_t gp = g_i; gp <= g_f; ++gp)
            if (gp < gs_i or gp > gs_f)
              b[i + g] += chi * xs->nu_sigma_f[gp] *
                          phi[uk_map + gp] * volume;
        }
        else
        {
          const double chi_p = xs->chi_prompt[g];
          for (size_t gp = gs_i; gp <= gs_f; ++gp)
            if (gp < gs_i or gp > gs_f)
              b[i + g] += chi_p * xs->nu_prompt_sigma_f[gp] *
                          phi[uk_map + gp] * volume;

          for (size_t j = 0; j < xs->n_precursors; ++j)
          {
            const double chi_d = xs->chi_delayed[g][j];
            const double gamma = xs->precursor_yield[j];
            for (size_t gp = gs_i; gp < gs_f; ++gp)
              if (gp < gs_i or gp > gs_f)
                b[i + g] += chi_d * gamma *
                            xs->nu_delayed_sigma_f[gp] *
                            phi[uk_map + gp] * volume;
          }
        }
      }

    }//for group

    //============================== Cross-group coupling
    if (solution_technique == SolutionTechnique::GROUPSET_WISE)
    {
      for (size_t g = gs_i; g <= gs_f; ++g)
      {
        //==================== Scattering term
        for (size_t gp = g_i; gp <= g_f; ++gp)
          b[i + g] += xs->transfer_matrices[0][g][gp] *
                      phi[uk_map + gp] * volume;

        //==================== Fission term
        if (xs->is_fissile)
        {
          //========== Total fission
          if (not use_precursors)
          {
            const double chi = xs->chi[g];
            for (size_t gp = g_i; gp <= g_f; ++gp)
              b[i + g] += chi * xs->nu_sigma_f[gp] *
                          phi[uk_map + gp] * volume;
          }

          //========== Prompt + delayed fission
          else
          {
            const double chi_p = xs->chi_prompt[g];
            for (size_t gp = g_i; gp <= g_f; ++gp)
              b[i + g] += chi_p * xs->nu_prompt_sigma_f[gp] *
                          phi[uk_map + gp] * volume;

            for (size_t j = 0; j < xs->n_precursors; ++j)
            {
              const double chi_d = xs->chi_delayed[g][j];
              const double gamma = xs->precursor_yield[j];
              for (size_t gp = g_i; gp <= g_f; ++gp)
                b[i + g] += chi_d * gamma *
                            xs->nu_delayed_sigma_f[gp] *
                            phi[uk_map + gp] * volume;
            }
          }
        }//if fissile
      }//for group
    }//if groupset-wise solver

    //================================================== Loop over faces
    for (const auto& face : cell->faces)
    {
      //============================== Boundary faces
      if (not face.neighbor_id)
      {
        const auto bndry_id = face.neighbor_id;
        const auto& bndry_type = boundary_info[bndry_id].first;

        //==================== Dirichlet boundary term
        if (bndry_type == BoundaryType::DIRICHLET)
        {
          const double d_pf = cell->centroid.distance(face.centroid);
          for (size_t g = gs_i; g <= gs_f; ++g)
          {
            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<DirichletBoundary>(bndry);
            const double D = xs->diffusion_coeff[g];
            b[i + g] += D/d_pf * bc->value * face.area;
          }
        }

        //==================== Neumann boundary term
        else if (bndry_type == BoundaryType::NEUMANN)
        {
          for (size_t g = gs_i; g <= gs_f; ++g)
          {
            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<NeumannBoundary>(bndry);
            b[i + g] += bc->value * face.area;
          }
        }

        //==================== Robin boundary term
        else if (bndry_type == BoundaryType::MARSHAK or
                 bndry_type == BoundaryType::ROBIN)
        {
          const double d_pf = cell->centroid.distance(face.centroid);
          for (size_t g = gs_i; g <= gs_f; ++g)
          {
            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<RobinBoundary>(bndry);
            const double D = xs->diffusion_coeff[g];
            b[i + g] += D/(bc->b*D + bc->a*d_pf) * bc->f * face.area;
          }
        }
      }//if boundary face
    }//for faces
  }//for cells
}
