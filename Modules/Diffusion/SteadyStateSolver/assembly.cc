#include "steadystate_solver.h"


//############################################################
void diffusion::SteadyStateSolver::assemble_matrix()
{
  switch (discretization->type)
  {
    case SpatialDiscretizationMethod::FINITE_VOLUME:
    {
      fv_assemble_matrix();
      break;
    }
    default:
      break;
  }
}


//############################################################
void diffusion::SteadyStateSolver::assemble_rhs_vector()
{
  switch (discretization->type)
  {
    case SpatialDiscretizationMethod::FINITE_VOLUME:
    {
      fv_assemble_rhs_vector();
      break;
    }
    default:
      break;
    }
  }


//############################################################
void diffusion::SteadyStateSolver::fv_assemble_matrix()
{
  system_matrix *= 0.0; // Zero out the matrix

  // Loop over cells
  for (const auto& cell : mesh->cells)
  {
    // Get cell information
    const auto& volume = cell->volume;
    const auto& xs = material_xs[matid_to_xs_map[cell->material_id]];

    // Loop over groups, add interaction term
    const int i = cell->id * n_groups; // starting index for cell
    for (int g = 0; g < n_groups; ++g)
      system_matrix[i + g][i + g] = xs->sigma_t[g] * volume;

    // Loop over faces
    for (const auto& face : cell->faces)
    {
      // Handle interior faces
      if (face.has_neighbor)
      {
        // Get neighbor information
        const auto& nbr_cell = mesh->cells[face.neighbor_id];
        const auto& nbr_xs = material_xs[matid_to_xs_map[nbr_cell->material_id]];

        // Get diffusion coefficients
        const auto& D = xs->diffusion_coeff;
        const auto& D_nbr = nbr_xs->diffusion_coeff;

        // Compute distances to neighbor and face centroid
        double d_n = cell->centroid.distance(nbr_cell->centroid);
        double d_f = cell->centroid.distance(face.centroid);

        // Weighting factor for harmonic mean
        double w = d_f / d_n;

        // Loop over groups
        const int j = nbr_cell->id * n_groups; // starting index for neighbor
        for (int g = 0; g < n_groups; ++g)
        {
          // Compute face diffusion coefficient
          double D_face = 1.0 / (w/D[g] + (1.0 - w)/D_nbr[g]);
          double value = D_face/d_n * face.area;
          system_matrix[i + g][i + g] += value;
          system_matrix[i + g][j + g] -= value;
        }
      }//if has neighbor

      // Apply boundary conditions
      else
      {
        double d_f = cell->centroid.distance(face.centroid);
        const int bndry_id = face.neighbor_id;
        auto bndry_type = boundary_info[bndry_id].first;

        // Dirichlet type boundaries
        if (bndry_type == BoundaryType::DIRICHLET or
            bndry_type == BoundaryType::ZERO_FLUX)
        {
          const auto& D = xs->diffusion_coeff;
          for (int g = 0; g < n_groups; ++g)
          {
            double value = D[g] / d_f * face.area;
            system_matrix[i + g][i + g] += value;
          }
        }//if Dirichlet boundary

        // Robin type boundaries
        else if (bndry_type == BoundaryType::ROBIN or
                 bndry_type == BoundaryType::VACUUM or
                 bndry_type == BoundaryType::MARSHAK)
        {
          const auto& D = xs->diffusion_coeff;
          for (int g = 0; g < n_groups; ++g)
          {
            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<RobinBoundary>(bndry);

            double value = bc->a*D[g] / (bc->a*d_f + bc->b*D[g]) * face.area;
            system_matrix[i + g][i + g] += value;
          }
        }//if Robin boundary
      }//if boundary face
    }//for face
  }//for cell

  // Contsturct the matrix
}


//############################################################
void diffusion::SteadyStateSolver::fv_assemble_rhs_vector()
{
  system_rhs *= 0.0; // Zero out the vector

  // Loop over cells
  for (const auto& cell : mesh->cells)
  {
    // Get source, if no source, go to next cell
    const int src_id = matid_to_src_map[cell->material_id];
    if (src_id < 0) continue;
    const auto& src = material_src[src_id];

    // Loop over groups
    const int i = cell->id * n_groups;
    for (int g = 0; g < n_groups; ++g)
      system_rhs[i + g] += src->values[g] * cell->volume;

    // Loop over faces to apply boundary values
    for (const auto& face : cell->faces)
    {
      if (not face.has_neighbor)
      {
        const int bndry_id = face.neighbor_id;
        const auto bndry_type = boundary_info[bndry_id].first;

        // Handle Neumann boundaries
        if (bndry_type == BoundaryType::NEUMANN)
        {
          for (int g = 0; g < n_groups; ++g)
          {
            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<NeumannBoundary>(bndry);

            system_rhs[i + g] += bc->value * face.area;
          }
        }//if Neumann

        // Handle Dirichlet and Robin boundaries
        else if (bndry_type == BoundaryType::DIRICHLET or
                 bndry_type == BoundaryType::ROBIN or
                 bndry_type == BoundaryType::MARSHAK)
        {
          // Get cell information
          double d_f = cell->centroid.distance(face.centroid);
          const auto& xs = material_xs[matid_to_xs_map[cell->material_id]];
          const auto& D = xs->diffusion_coeff;

          // Dirichlet boundaries
          if (bndry_type == BoundaryType::DIRICHLET)
          {
            for (int g = 0; g < n_groups; ++g)
            {
              const auto& bndry = boundaries[bndry_id][g];
              const auto bc =
                  std::static_pointer_cast<DirichletBoundary>(bndry);

              system_rhs[i + g] += D[g]/d_f * bc->value * face.area;
            }
          }//if Dirichlet

          // Robin boundaries
          else
          {
            for (int g = 0; g < n_groups; ++g)
            {
              const auto& bndry = boundaries[bndry_id][g];
              const auto bc = std::static_pointer_cast<RobinBoundary>(bndry);

              system_rhs[i + g] +=
                  D[g]/(bc->b*D[g] + bc->a*d_f) * bc->f * face.area;
            }
          }//if Robin
        }//if Dirichlet or Robin
      }//if boundary face
    }//for face
  }//for cell
}
