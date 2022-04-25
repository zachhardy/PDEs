#include "steadystate_solver.h"

/**
 * \brief Assemble the full multigroup diffusion operator using the finite
 *        volume method.
 *
 * This routine first clears the \p system_matrix. It then assembles the full
 * multigroup operator, including the group-wise total interaction, diffusion,
 * and boundary terms as well as the cross-group scattering and fission term.
 */
void neutron_diffusion::SteadyStateSolver::assemble_fv_matrix()
{
  system_matrix *= 0.0; // Clear the matrix

  // Build the system cell by cell
  for (const auto& cell : mesh->cells)
  {
    double volume = cell->volume;
    const auto& xs = material_xs[matid_to_xs_map[cell->material_id]];
    const size_t i = cell->id * n_groups;

    // Assemble group-wise cell terms
    for (size_t g = 0; g < n_groups; ++g)
    {
      // Total interaction term
      system_matrix[i + g][i + g] += xs->sigma_t[g] * volume;

      // Assemble scattering term
      for (size_t gp = 0; gp < n_groups; ++gp)
        system_matrix[i + g][i + gp] -=
            xs->transfer_matrices[0][g][gp] * volume;

      // Assemble fission term
      if (xs->is_fissile)
      {
        double chi = (not use_precursors)? xs->chi[g] : xs->chi_prompt[g];

        for (size_t gp = 0; g < n_groups; ++gp)
        {
          double nu_sig_f = (not use_precursors)? xs->nu_sigma_f[gp] :
                                                  xs->nu_prompt_sigma_f[gp];
          system_matrix[i + g][i + gp] -= chi * nu_sig_f * volume;
        }
      }//if fissile
    }//for group

    // Assemble group-wise face terms
    for (const auto& face : cell->faces)
    {
      // Diffusion across interior faces
      if (face.has_neighbor)
      {
        const auto& nbr_cell = mesh->cells[face.neighbor_id];
        const size_t j = nbr_cell->id * n_groups;

        const int nbr_xs_id = matid_to_xs_map[nbr_cell->material_id];
        const auto& nbr_xs = material_xs[nbr_xs_id];

        // Geometric quantities
        double d_pn = cell->centroid.distance(nbr_cell->centroid);
        double d_pf = cell->centroid.distance(face.centroid);
        double w = d_pf / d_pn; // Weighting factor for harmonic means

        // Diffusion coefficients
        const auto& D = xs->diffusion_coeff;
        const auto& D_nbr = nbr_xs->diffusion_coeff;

        for (size_t g = 0; g < n_groups; ++g)
        {
          // Effective diffusion coefficient on face
          double D_eff = 1.0 / (w/D[g] + (1.0 - w)/D_nbr[g]);

          system_matrix[i + g][i + g] += D_eff/d_pn * face.area;
          system_matrix[i + g][j + g] -= D_eff/d_pn * face.area;
        }//for group
      }//if interior face

      // Assemble group-wise boundary terms
      else
      {
        const auto bndry_id = face.neighbor_id;
        const auto bndry_type = boundary_info[bndry_id].first;

        // Handle Dirichlet boundaries
        if (bndry_type == BoundaryType::ZERO_FLUX or
            bndry_type == BoundaryType::DIRICHLET)
        {
          double d_pf = cell->centroid.distance(face.centroid);
          const auto& D = xs->diffusion_coeff;

          for (size_t g = 0; g < n_groups; ++g)
            system_matrix[i + g][i + g] += D[g]/d_pf * face.area;
        }

        // Handle Robin boundaries
        if (bndry_type == BoundaryType::VACUUM or
            bndry_type == BoundaryType::MARSHAK or
            bndry_type == BoundaryType::ROBIN)
        {
          double d_pf = cell->centroid.distance(face.centroid);
          const auto& D = xs->diffusion_coeff;

          for (size_t g = 0; g < n_groups; ++g)
          {
            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<RobinBoundary>(bndry);
            system_matrix[i + g][i + g] +=
                bc->a*D[g]/(bc->b*D[g] + bc->a*d_pf) * face.area;
          }
        }
      }//if boundary face
    }//for face
  }//for cell
}

//######################################################################

/**
 * \brief Assemble the full multigroup right-hand side vector using the
 *        finite volume method.
 *
 * This routine first clears the \p system_rhs. It then assembles the full
 * multigroup right-hand side vector, including the multigroup inhomogeneous
 * source and boundary sources.
 */
void neutron_diffusion::SteadyStateSolver::assemble_fv_rhs()
{
  system_rhs *= 0.0; // Clear the rhs

  // Assemble cell by cell
  for (const auto& cell : mesh->cells)
  {
    double volume = cell->volume;
    const size_t i = cell->id * n_groups;

    const int src_id = matid_to_src_map[cell->material_id];
    if (src_id < 0)
      continue;
    const auto src = material_src[src_id];

    // Contribute group-wise inhomogeneous source
    for (size_t g = 0; g < n_groups; ++g)
      system_rhs[i + g] += src->values[g] * volume;

    // Assemble group-wise boundary terms
    for (const auto& face : cell->faces)
    {
      // Only consider boundary faces
      if (not face.has_neighbor)
      {
        const size_t bndry_id = face.neighbor_id;
        const auto bndry_type = boundary_info[bndry_id].first;

        // Handle Dirichlet boundaries
        if (bndry_type == BoundaryType::DIRICHLET)
        {
          double d_pf = cell->centroid.distance(face.centroid);

          const auto& xs = material_xs[matid_to_xs_map[cell->material_id]];
          const auto& D = xs->diffusion_coeff;

          for (size_t g = 0; g < n_groups; ++g)
          {
            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<DirichletBoundary>(bndry);
            system_rhs[i + g] += D[g]/d_pf * bc->value * face.area;
          }
        }//if Dirichlet

        // Handle Neumann boundaries
        else if (bndry_type == BoundaryType::NEUMANN)
        {
          for (size_t g = 0; g < n_groups; ++g)
          {
            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<NeumannBoundary>(bndry);
            system_rhs[i + g] += bc->value * face.area;
          }
        }//if Neumann

        // Handle Robin boundaries
        else if (bndry_type == BoundaryType::MARSHAK or
                 bndry_type == BoundaryType::ROBIN)
        {
          double d_pf = cell->centroid.distance(face.centroid);

          const auto& xs = material_xs[matid_to_xs_map[cell->material_id]];
          const auto& D = xs->diffusion_coeff;

          for (size_t g = 0; g < n_groups; ++g)
          {
            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<RobinBoundary>(bndry);
            system_rhs[i + g] +=
                D[g]/(bc->b*D[g] + bc->a*d_pf) * bc->f * face.area;
          }
        }//if Robin
      }
    }//for face
  }//for cell
}
