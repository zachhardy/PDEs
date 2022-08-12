#include "steadystate_solver.h"


using namespace Math;
using namespace NeutronDiffusion;


void
SteadyStateSolver::
assemble_matrix(AssemblerFlags assembler_flags)
{
  const bool assemble_scatter = (assembler_flags & ASSEMBLE_SCATTER);
  const bool assemble_fission = (assembler_flags & ASSEMBLE_FISSION);

  A = 0.0;

  // Loop over cells
  for (const auto& cell : mesh->cells)
  {
    const auto& volume = cell.volume;
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];
    const auto i = n_groups * cell.id;

    const auto* sig_t = cellwise_xs[cell.id].sigma_t.data();
    const auto* D = xs->diffusion_coeff.data();
    const auto* B = xs->buckling.data();

    // Loop over groups
    for (unsigned int gr = 0; gr < n_groups; ++gr)
    {
      const auto g = groups[gr];

      //========================================
      // Total interaction term + buckling
      //========================================

      A.add(i + gr, i + gr, (sig_t[g] + D[g]*B[g]) * volume);

      //========================================
      // Scattering term
      //========================================

      if (assemble_scatter)
      {
        const auto* sig_s = xs->transfer_matrices[0][g].data();
        for (unsigned int gpr = 0; gpr < n_groups; ++gpr)
          A.add(i + gr, i + gpr, -sig_s[groups[gpr]] * volume);
      }

      //========================================
      // Fission term
      //========================================

      if (xs->is_fissile && assemble_fission)
      {
        // Total fission
        if (not use_precursors)
        {
          const auto chi = xs->chi[g];
          const auto* nu_sigf = xs->nu_sigma_f.data();

          for (unsigned int gpr = 0; gpr < n_groups; ++gpr)
            A.add(i + gr, i + gpr, -chi * nu_sigf[groups[gpr]] * volume);
        }

        // Prompt + delayed fission
        else
        {
          const auto chi_p = xs->chi_prompt[g];
          const auto* chi_d = xs->chi_delayed[g].data();
          const auto* nup_sigf = xs->nu_prompt_sigma_f.data();
          const auto* nud_sigf = xs->nu_delayed_sigma_f.data();
          const auto* gamma = xs->precursor_yield.data();

          for (unsigned int gpr = 0; gpr < n_groups; ++gpr)
          {
            double f = chi_p*nup_sigf[groups[gpr]];
            for (unsigned int j = 0; j < xs->n_precursors; ++j)
              f += chi_d[j] * gamma[j] * nud_sigf[groups[gpr]];
            A.add(i + gr, i + gpr, -f * volume);
          }
        }
      }//if fissile
    }//for group

    // Loop over faces
    for (const auto& face : cell.faces)
    {
      //========================================
      // Diffusion term on interior faces
      //========================================

      if (face.has_neighbor)
      {
        // Get neighbor info
        const auto& nbr_cell = mesh->cells[face.neighbor_id];
        const auto nbr_xs_id = matid_to_xs_map[nbr_cell.material_id];
        const auto& nbr_xs = material_xs[nbr_xs_id];
        const auto* D_nbr = nbr_xs->diffusion_coeff.data();

        const auto j = n_groups * nbr_cell.id;

        // Geometric quantities
        const auto d_pf = cell.centroid.distance(face.centroid);
        const auto d_pn = cell.centroid.distance(nbr_cell.centroid);
        const auto w = d_pf/d_pn; // harmonic mean weighting factor

        for (unsigned int gr = 0; gr < n_groups; ++gr)
        {
          const auto g = groups[gr];

          const double D_eff = 1.0/(w/D[g] + (1.0 - w)/D_nbr[g]);
          const double value = D_eff/d_pn * face.area;

          A.add(i + gr, i + gr, value);
          A.add(i + gr, j + gr, -value);
        }
      }//if interior face

      //========================================
      // Boundary terms
      //========================================

      else
      {
        const auto bndry_id = face.neighbor_id;
        const auto bndry_type = boundary_info[bndry_id].first;

        //========================================
        // Dirichlet boundary term
        //========================================

        if (bndry_type == BoundaryType::ZERO_FLUX or
            bndry_type == BoundaryType::DIRICHLET)
        {
          const auto d_pf = cell.centroid.distance(face.centroid);
          for (unsigned int gr = 0; gr < n_groups; ++gr)
            A.add(i + gr, i + gr, D[groups[gr]]/d_pf * face.area);
        }

        //========================================
        // Robin boundary term
        //========================================

        else if (bndry_type == BoundaryType::VACUUM or
                 bndry_type == BoundaryType::MARSHAK or
                 bndry_type == BoundaryType::ROBIN)
        {
          const auto d_pf = cell.centroid.distance(face.centroid);
          for (unsigned int gr = 0; gr < n_groups; ++gr)
          {
            const auto g = groups[gr];

            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<RobinBoundary>(bndry);

            const auto value = bc->a*D[g]/(bc->b*D[g] + bc->a*d_pf);
            A.add(i + gr, i + gr, value * face.area);
          }
        }
      }//if boundary face
    }//for face
  }//for cell
}
