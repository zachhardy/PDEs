#include "steadystate_solver_fv.h"


using namespace pdes;
using namespace Math;

/**
* Assemble the matrix for the specified \p groupset.
 *
 * If solving the full system, this routine assembles the full multigroup
 * operator for the groupset, including off-diagonal scattering and fission
 * coupling terms. Otherwise, this routine assembles the within-group system
 * for the groups within the groupset.
*/
void NeutronDiffusion::SteadyStateSolver_FV::
assemble_matrix(Groupset& groupset)
{
  Matrix A(groupset.matrix.n_rows(),
           groupset.matrix.n_cols(), 0.0);

  // Get groupset range
  const auto n_gsg = groupset.groups.size();
  const auto gs_i = groupset.groups.front();
  const auto gs_f = groupset.groups.back();

  //================================================== Loop over cells
  for (const auto& cell : mesh->cells)
  {
    const double volume = cell.volume;
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];
    const size_t i = cell.id * n_gsg;

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
        {
          double value = xs->transfer_matrices[0][g][gp] * volume;
          A[i + g][i + gp] -= value;
        }

        //==================== Fission term
        if (xs->is_fissile)
        {
          //========== Total fission
          if (not use_precursors)
          {
            const double chi = xs->chi[g];
            for (size_t gp = gs_i; gp <= gs_f; ++gp)
            {
              double value = chi * xs->nu_sigma_f[gp] * volume;
              A[i + g][i + gp] -= value;
            }
          }

          //========== Prompt + delayed fission
          else
          {
            const double chi_p = xs->chi_prompt[g];
            for (size_t gp = gs_i; gp <= gs_f; ++gp)
            {
              double value = chi_p * xs->nu_prompt_sigma_f[gp] * volume;
              A[i + g][i + gp] -= value;
            }

            for (size_t j = 0; j < xs->n_precursors; ++j)
            {
              const double chi_d = xs->chi_delayed[g][j];
              const double gamma = xs->precursor_yield[j];
              for (size_t gp = gs_i; gp <= gs_f; ++gp)
              {
                double value =
                    chi_d * gamma * xs->nu_delayed_sigma_f[gp] * volume;
                A[i + g][i + gp] -= value;
              }
            }//for precursors
          }
        }//if fissile
      }//for group
    }//if full system

    //================================================== Loop over faces
    for (const auto& face : cell.faces)
    {
      //============================== Interior faces
      if (face.has_neighbor)
      {
        // Get neighbor info
        const auto& nbr_cell = mesh->cells[face.neighbor_id];
        const size_t j = nbr_cell.id * n_gsg;

        const int nbr_xs_id = matid_to_xs_map[nbr_cell.material_id];
        const auto& nbr_xs = material_xs[nbr_xs_id];

        // Geometric quantities
        const double d_pf = cell.centroid.distance(face.centroid);
        const double d_pn = cell.centroid.distance(nbr_cell.centroid);
        const double w = d_pf / d_pn; // harmonic mean weighting factor

        //==================== Diffusion term
        for (size_t g = gs_i; g <= gs_f; ++g)
        {
          const double D = xs->diffusion_coeff[g];
          const double D_nbr = nbr_xs->diffusion_coeff[g];
          const double D_eff = 1.0 / (w/D + (1.0 - w)/D_nbr);

          double value = D_eff / d_pn * face.area;
          A[i + g][i + g] += value;
          A[i + g][j + g] -= value;
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
          const double d_pf = cell.centroid.distance(face.centroid);
          for (size_t g = gs_i; g <= gs_f; ++g)
          {
            const double D = xs->diffusion_coeff[g];
            double value = D / d_pf * face.area;
            A[i + g][i + g] += value;
          }
        }

        //==================== Robin boundary term
        else if (bndry_type == BoundaryType::VACUUM or
                 bndry_type == BoundaryType::MARSHAK or
                 bndry_type == BoundaryType::ROBIN)
        {
          const double d_pf = cell.centroid.distance(face.centroid);
          for (size_t g = gs_i; g <= gs_f; ++g)
          {
            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<RobinBoundary>(bndry);
            const double D = xs->diffusion_coeff[g];
            double value = bc->a*D/(bc->b*D + bc->a*d_pf) * face.area;
            A[i + g][i + g] += value;
          }
        }
      }//if boundary face
    }//for face
  }//for cell

  groupset.matrix = A;
}


