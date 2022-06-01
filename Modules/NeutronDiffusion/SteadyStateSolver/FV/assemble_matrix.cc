#include "../steadystate_solver.h"


using namespace Math;
using namespace NeutronDiffusion;


void
SteadyStateSolver::
fv_assemble_matrix(Groupset& groupset, AssemblerFlags assembler_flags)
{
  const bool assemble_scatter = (assembler_flags & ASSEMBLE_SCATTER);
  const bool assemble_fission = (assembler_flags & ASSEMBLE_FISSION);

  SparseMatrix& A = groupset.A = 0.0;

  // Get groupset range
  const size_t n_gsg = groupset.groups.size();

  //================================================== Loop over cells
  for (const auto& cell : mesh->cells)
  {
    const double volume = cell.volume;
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];
    const size_t i = cell.id * n_gsg;

    //============================== Loop over groups
    for (size_t gr = 0; gr < n_gsg; ++gr)
    {
      const size_t ig = i + gr;
      const size_t g = groupset.groups[gr];

      //==================== Total interaction term
      A.add(ig, ig, xs->sigma_t[g] * volume);

      if (assemble_scatter)
      {
        //==================== Scattering term
        const double* sig_s = &xs->transfer_matrices[0][g][0];

        for (size_t gpr = 0; gpr < n_gsg; ++gpr)
        {
          const size_t igp = i + gpr;
          const size_t gp = groupset.groups[gpr];
          A.add(ig, igp, -sig_s[gp] * volume);
        }
      }


      //==================== Fission term
      if (xs->is_fissile && assemble_fission)
      {
        //========== Total fission
        if (not use_precursors)
        {
          const double chi = xs->chi[g];
          const double* nu_sigf = &xs->nu_sigma_f[0];

          for (size_t gpr = 0; gpr < n_gsg; ++gpr)
          {
            const size_t igp = i + gpr;
            const size_t gp = groupset.groups[gpr];
            A.add(ig, igp, -chi * nu_sigf[gp] * volume);
          }
        }

          //========== Prompt + delayed fission
        else
        {
          const double chi_p = xs->chi_prompt[g];
          const double* chi_d = &xs->chi_delayed[g][0];
          const double* nup_sigf = &xs->nu_prompt_sigma_f[0];
          const double* nud_sigf = &xs->nu_delayed_sigma_f[0];
          const double* gamma = &xs->precursor_yield[0];

          for (size_t gpr = 0; gpr < n_gsg; ++gpr)
          {
            const size_t igp = i + gpr;
            const size_t gp = groupset.groups[gpr];

            double f = chi_p * nup_sigf[gp];
            for (size_t j = 0; j < xs->n_precursors; ++j)
              f += chi_d[j] * gamma[j] * nud_sigf[gp];
            A.add(ig, igp, -f * volume);
          }
        }
      }//if fissile
    }//for group

    //============================================= Loop over faces
    for (const auto& face : cell.faces)
    {
      //============================== Interior faces
      if (face.has_neighbor)
      {
        // Get neighbor info
        const auto& nbr_cell = mesh->cells[face.neighbor_id];
        const int nbr_xs_id = matid_to_xs_map[nbr_cell.material_id];
        const auto& nbr_xs = material_xs[nbr_xs_id];

        const size_t j = nbr_cell.id * n_gsg;

        // Geometric quantities
        const double d_pf = cell.centroid.distance(face.centroid);
        const double d_pn = cell.centroid.distance(nbr_cell.centroid);
        const double w = d_pf / d_pn; // harmonic mean weighting factor

        //==================== Diffusion term
        const double* D = &xs->diffusion_coeff[0];
        const double* D_nbr = &nbr_xs->diffusion_coeff[0];

        for (size_t gr = 0; gr < n_gsg; ++gr)
        {
          const size_t ig = i + gr;
          const size_t jg = j + gr;
          const size_t g = groupset.groups[gr];

          const double D_eff = 1.0 / (w / D[g] + (1.0 - w) / D_nbr[g]);
          const double value = D_eff / d_pn * face.area;

          A.add(ig, ig, value);
          A.add(ig, jg, -value);
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
          const double* D = &xs->diffusion_coeff[0];
          const double d_pf = cell.centroid.distance(face.centroid);

          for (size_t gr = 0; gr < n_gsg; ++gr)
          {
            const size_t ig = i + gr;
            const size_t g = groupset.groups[gr];
            A.add(ig, ig, D[g] / d_pf * face.area);
          }
        }

          //==================== Robin boundary term
        else if (bndry_type == BoundaryType::VACUUM or
                 bndry_type == BoundaryType::MARSHAK or
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

            double value = bc->a * D[g] / (bc->b * D[g] + bc->a * d_pf);
            A.add(ig, ig, value * face.area);
          }
        }
      }//if boundary face
    }//for face
  }//for cell
}
