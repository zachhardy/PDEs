#include "steadystate_solver.h"

void neutron_diffusion::SteadyStateSolver::assemble_fv_matrix()
{
  math::Matrix& A = system_matrix;
  A *= 0.0;

  //================================================== Loop over cells
  for (const auto& cell : mesh->cells)
  {
    double volume = cell->volume;
    const auto& xs = material_xs[matid_to_xs_map[cell->material_id]];
    size_t i = cell->id * n_groups;

    // Total interaction term
    for (size_t g = 0; g < n_groups; ++g)
      A[i + g][i + g] += xs->sigma_t[g] * volume;

    //================================================== Group transfers
    if (solution_method == SolutionMethod::DIRECT)
    {
      for (size_t g = 0; g < n_groups; ++g)
      {
        // Scattering term
        for (size_t gp = 0; gp < n_groups; ++gp)
          A[i + g][i + gp] -= xs->transfer_matrices[0][g][gp] * volume;

        // Fission term
        if (xs->is_fissile)
        {
          // Total fission term
          if (not use_precursors)
          {
            double chi = xs->chi[g];
            for (size_t gp = 0; gp < n_groups; ++gp)
              A[i + g][i + gp] -= chi * xs->nu_sigma_f[gp] * volume;
          }

          // Prompt + delayed term
          else
          {
            double chi_p = xs->chi_prompt[g];
            for (size_t gp = 0; gp < n_groups; ++gp)
              A[i + g][i + gp] -= chi_p * xs->nu_prompt_sigma_f[gp] * volume;

            for (size_t j = 0; j < xs->n_precursors; ++j)
            {
              double chi_d = xs->chi_delayed[g][j];
              double gamma = xs->precursor_yield[j];
              for (size_t gp = 0; gp < n_groups; ++gp)
                A[i + g][i + gp] -=
                    chi_d * gamma * xs->nu_delayed_sigma_f[gp] * volume;
            }
          }
        }//if fissile
      }//for group
    }//if off-diagonal terms

    //================================================== Loop over faces
    for (const auto& face : cell->faces)
    {
      //================================================== Interior faces
      // Diffusion term
      if (face.has_neighbor)
      {
        // Get neighbor cell
        const auto& nbr_cell = mesh->cells[face.neighbor_id];
        int nbr_mat_id = nbr_cell->material_id;
        const auto& nbr_xs = material_xs[nbr_mat_id];
        size_t j = nbr_cell->id * n_groups;

        // Geometric quantities
        double d_pf = cell->centroid.distance(face.centroid);
        double d_pn = cell->centroid.distance(nbr_cell->centroid);
        double w = d_pf / d_pn; // weighting factor for harmonic mean

        for (size_t g = 0; g < n_groups; ++g)
        {
          // Diffusion coefficients
          double D = xs->diffusion_coeff[g];
          double D_nbr = nbr_xs->diffusion_coeff[g];
          double D_eff = 1.0 / (w/D + (1.0 - w)/D_nbr);

          A[i + g][i + g] += D_eff / d_pn * face.area;
          A[i + g][j +g] -= D_eff / d_pn * face.area;
        }
      }//if interior face

      //================================================== Boundary faces
      else
      {
        const size_t bndry_id = face.neighbor_id;
        const auto& bndry_type = boundary_info[bndry_id].first;

        // Dirichlet boundary term
        if (bndry_type == BoundaryType::ZERO_FLUX or
            bndry_type == BoundaryType::DIRICHLET)
        {
          double d_pf = cell->centroid.distance(face.centroid);
          for (size_t g = 0; g < n_groups; ++g)
          {
            double D = xs->diffusion_coeff[g];
            A[i + g][i + g] += D / d_pf * face.area;
          }
        }

        // Robin boundary term
        else if (bndry_type == BoundaryType::VACUUM or
                 bndry_type == BoundaryType::MARSHAK or
                 bndry_type == BoundaryType::ROBIN)
        {
          double d_pf = cell->centroid.distance(face.centroid);
          for (size_t g = 0; g < n_groups; ++g)
          {
            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<RobinBoundary>(bndry);

            double D = xs->diffusion_coeff[g];
            A[i + g][i + g] += bc->a*D/(bc->b*D + bc->a*d_pf) * face.area;
          }
        }
      }//if boundary face
    }//for face
  }//for cell
}

//######################################################################

void neutron_diffusion::SteadyStateSolver::set_fv_source()
{
  math::Vector& b = system_rhs;
  b *= 0.0;

  //================================================== Loop over cells
  for (const auto& cell : mesh->cells)
  {
    double volume = cell->volume;
    const auto& xs = material_xs[matid_to_xs_map[cell->material_id]];
    size_t i = cell->id * n_groups;

    // Inhomogeneous source
    int src_id = matid_to_src_map[cell->material_id];
    if (src_id >= 0)
    {
      const auto& src = material_src[src_id];
      for (size_t g = 0; g < n_groups; ++g)
        b[i + g] += src->values[g] * volume;
    }

    // ================================================== Group transfers
    if (solution_method == SolutionMethod::ITERATIVE)
    {
      for (size_t g = 0; g < n_groups; ++g)
      {
        // Scattering term
        for (size_t gp = 0; gp < n_groups; ++gp)
          b[i + g] += xs->transfer_matrices[0][g][gp] *
                    phi[i + gp] * volume;

        // Fission term
        if (xs->is_fissile)
        {
          // Total fission
          if (not use_precursors)
          {
            double chi = xs->chi[g];
            for (size_t gp = 0; gp < n_groups; ++gp)
              b[i + g] += chi * xs->nu_sigma_f[gp] *
                        phi[i + gp] * volume;

            for (size_t j = 0; j < xs->n_precursors; ++j)
            {
              double chi_d = xs->chi_delayed[g][j];
              double gamma = xs->precursor_yield[j];
              for (size_t gp = 0; gp < n_groups; ++gp)
                b[i + g] += chi_d * gamma *
                          xs->nu_delayed_sigma_f[gp] *
                          phi[i + gp] * volume;
            }
          }
        }//if fissile
      }//for group
    }//if including group transfers

    //================================================== Loop over faces
    for (const auto& face : cell->faces)
    {
      //================================================== Boundary faces
      if (not face.has_neighbor)
      {
        const auto bndry_id = face.neighbor_id;
        const auto bndry_type = boundary_info[bndry_id].first;

        // Dirichlet boundary term
        if (bndry_type == BoundaryType::DIRICHLET)
        {
          double d_pf = cell->centroid.distance(face.centroid);
          for (size_t g = 0; g < n_groups; ++g)
          {
            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<DirichletBoundary>(bndry);

            double D = xs->diffusion_coeff[g];
            b[i + g] += D / d_pf * bc->value * face.area;
          }
        }

        // Neumann boundary term
        else if (bndry_type == BoundaryType::NEUMANN)
        {
          for (size_t g = 0; g < n_groups; ++g)
          {
            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<NeumannBoundary>(bndry);
            b[i + g] += bc->value * face.area;
          }
        }

        // Robin boundary term
        else if (bndry_type == BoundaryType::MARSHAK or
                 bndry_type == BoundaryType::ROBIN)
        {
          double d_pf = cell->centroid.distance(face.centroid);
          for (size_t g = 0; g < n_groups; ++g)
          {
            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<RobinBoundary>(bndry);

            double D = xs->diffusion_coeff[g];
            b[i + g] += D/(bc->b*D + bc->a*d_pf) * bc->f * face.area;
          }
        }
      }//if boundary face
    }//for face
  }//for cell
}
