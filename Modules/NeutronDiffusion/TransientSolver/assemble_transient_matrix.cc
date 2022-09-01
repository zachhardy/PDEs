#include "transient_solver.h"

#include <iomanip>

using namespace NeutronDiffusion;


void
TransientSolver::
assemble_transient_matrix(AssemblerFlags assembler_flags)
{
  const bool assemble_scatter = (assembler_flags & ASSEMBLE_SCATTER);
  const bool assemble_fission = (assembler_flags & ASSEMBLE_FISSION);

  A = 0.0;

  // Get effective time step size
  const auto eff_dt = effective_time_step();

  // Loop over cells
  for (const auto& cell : mesh->cells)
  {
    const auto volume = cell.volume;
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];
    const auto i = n_groups * cell.id ;

    const auto* sig_t = cellwise_xs[cell.id].sigma_t.data();
    const auto* D = xs->diffusion_coeff.data();
    const auto* B = xs->buckling.data();
    const auto* inv_vel = xs->inv_velocity.data();

    // Loop over groups
    for (unsigned int g = 0; g < n_groups; ++g)
    {
      //=======================================================
      // Total interaction + buckling + time derivative term
      //=======================================================

      double entry = 0.0;
      entry += sig_t[g]; // total interaction
      entry += D[g] * B[g]; // buckling
      entry += inv_vel[g]/eff_dt; // time-derivative
      A.add(i + g, i + g, entry * volume);

      //========================================
      // Scattering term
      //========================================

      if (assemble_scatter)
      {
        const auto* sig_s = xs->transfer_matrices[0][g].data();
        for (unsigned int gp = 0; gp < n_groups; ++gp)
          A.add(i + g, i + gp, -sig_s[gp] * volume);
      }//if scattering

      //========================================
      // Fission term
      //========================================

      if (assemble_fission && xs->is_fissile)
      {
        // Total fission
        if (not use_precursors)
        {
          const auto chi = xs->chi[g];
          const auto* nu_sigf = xs->nu_sigma_f.data();

          for (unsigned int gp = 0; gp < n_groups; ++gp)
            A.add(i + g, i + gp, -chi * nu_sigf[gp] * volume);
        }//if no precursors

        //========== Prompt + delayed fission
        else
        {
          //===== Prompt
          const auto chi_p = xs->chi_prompt[g];
          const auto* nup_sigf = xs->nu_prompt_sigma_f.data();
          for (unsigned int gp = 0; gp < n_groups; ++gp)
            A.add(i + g, i + gp, -chi_p * nup_sigf[gp] * volume);

          //===== Delayed
          if (not lag_precursors)
          {
            const auto* chi_d = xs->chi_delayed[g].data();
            const auto* nud_sigf = xs->nu_delayed_sigma_f.data();
            const auto* lambda = xs->precursor_lambda.data();
            const auto* gamma = xs->precursor_yield.data();

            // This section of code computes the matrix coefficient
            // corresponding to the precursor substitution term arising when
            // precursors are treated implicitly. This term is similar to the
            // normal fission term with the exception of having an inner sum
            // over all precursor species. The first loop computes sum over
            // precursor species. This provides a factor which is multiplied
            // by nu_delayed_sigma_f when contributing to the matrix. The
            // coefficient is computed ahead of time for the sake of efficiency.
            // Using nested loops would require either this factor to be
            // computed \p n_gsg times or for many more expensive sparse matrix
            // add operations. */
            double coeff= 0.0;
            for (unsigned int j = 0; j < xs->n_precursors; ++j)
              coeff += chi_d[j] * lambda[j] * gamma[j] * eff_dt /
                       (1.0 + eff_dt*lambda[j]);

            // Contribute the delayed fission term to the matrix.
            for (unsigned int gp = 0; gp < n_groups; ++gp)
              A.add(i + g, i + gp, -coeff * nud_sigf[gp] * volume);
          }//if not lag precursors
        }//if prompt + delayed fission
      }//if fissile
    }//for group

    for (const auto& face : cell.faces)
    {
      //========================================
      // Diffusion term
      //========================================

      if (face.has_neighbor)
      {
        // Get neighbor info
        const auto& nbr_cell = mesh->cells[face.neighbor_id];
        const auto nbr_xs_id = matid_to_xs_map[nbr_cell.material_id];
        const auto& nbr_xs = material_xs[nbr_xs_id];
        const auto* D_nbr = nbr_xs->diffusion_coeff.data();

        const auto j = n_groups * nbr_cell.id;

        // Geometric data
        const auto d_pf = cell.centroid.distance(face.centroid);
        const auto d_pn = cell.centroid.distance(nbr_cell.centroid);
        const auto w = d_pf / d_pn; // harmonic mean weight factor

        for (unsigned int g = 0; g < n_groups; ++g)
        {
          const auto D_eff = 1.0/(w/D[g] + (1.0 - w)/D_nbr[g]);

          A.add(i + g, i + g, D_eff / d_pn * face.area);
          A.add(i + g, j + g, -D_eff / d_pn * face.area);
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

        if (bndry_type == BoundaryType::ZERO_FLUX ||
            bndry_type == BoundaryType::DIRICHLET)
        {
          const auto d_pf = cell.centroid.distance(face.centroid);
          for (unsigned int g = 0; g < n_groups; ++g)
            A.add(i + g, i + g, D[g] / d_pf * face.area);
        }//if Dirichlet

        //========================================
        // Robin boundary term
        //========================================

        else if (bndry_type == BoundaryType::VACUUM ||
                 bndry_type == BoundaryType::MARSHAK ||
                 bndry_type == BoundaryType::ROBIN)
        {
          const auto d_pf = cell.centroid.distance(face.centroid);
          for (unsigned int g = 0; g < n_groups; ++g)
          {
            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<RobinBoundary>(bndry);

            double val = bc->a*D[g]/(bc->b*D[g] + bc->a*d_pf);
            A.add(i + g, i + g, val * face.area);
          }//for g
        }//if Robin
      }//if boundary face
    }//for face
  }//for cell
}


void
TransientSolver::
rebuild_matrix()
{
  if (algorithm == Algorithm::DIRECT)
    assemble_transient_matrix(ASSEMBLE_SCATTER | ASSEMBLE_FISSION);
  else
    assemble_transient_matrix(NO_ASSEMBLER_FLAGS);
  linear_solver->set_matrix(A);
}
