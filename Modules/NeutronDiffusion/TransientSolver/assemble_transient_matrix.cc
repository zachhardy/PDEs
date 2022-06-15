#include "transient_solver.h"


using namespace NeutronDiffusion;


void
TransientSolver::assemble_transient_matrix(Groupset& groupset,
                                           AssemblerFlags assembler_flags)
{
  const bool assemble_scatter = (assembler_flags & ASSEMBLE_SCATTER);
  const bool assemble_fission = (assembler_flags & ASSEMBLE_FISSION);

  SparseMatrix& A = groupset.A = 0.0;

  // Get groupset range
  const size_t n_gsg = groupset.groups.size();

  // Get effective time step size
  const double eff_dt = effective_time_step();

  //======================================== Loop over cells
  for (const auto& cell : mesh->cells)
  {
    const double volume = cell.volume;
    const auto& xs = material_xs[matid_to_xs_map[cell.material_id]];
    const size_t uk_map = cell.id * n_gsg;

    const double* sig_t = xs->sigma_t.data();
    const double* inv_vel = xs->inv_velocity.data();

    //============================== Loop over groups
    for (size_t gr = 0; gr < n_gsg; ++gr)
    {
      const size_t g = groupset.groups[gr];

      //==================== Within group terms
      double entry = 0.0;
      entry += sig_t[g]; // total interaction
      entry += inv_vel[g]/eff_dt; // time-derivative
      A.add(uk_map + gr, uk_map + gr, entry * volume);

      //=============== Scattering term
      if (assemble_scatter)
      {
        const double* sig_s = xs->transfer_matrices[0][g].data();

        for (size_t gpr = 0; gpr < n_gsg; ++gpr)
        {
          const size_t gp = groupset.groups[gpr];
          A.add(uk_map + gr, uk_map + gpr, -sig_s[gp] * volume);
        }
      }//if scattering

      //=============== Fission term
      if (assemble_fission && xs->is_fissile)
      {
        //========== Total fission
        if (not use_precursors)
        {
          const double chi = xs->chi[g];
          const double* nu_sigf = xs->nu_sigma_f.data();

          for (size_t gpr = 0; gpr < n_gsg; ++gpr)
          {
            const size_t gp = groupset.groups[gpr];
            A.add(uk_map + gr, uk_map + gpr, -chi*nu_sigf[gp] * volume);
          }
        }//if no precursors

        //========== Prompt + delayed fission
        else
        {
          //===== Prompt
          const double chi_p = xs->chi_prompt[g];
          const double* nup_sigf = xs->nu_prompt_sigma_f.data();

          for (size_t gpr = 0; gpr < n_gsg; ++gpr)
          {
            const size_t gp = groupset.groups[gpr];
            A.add(uk_map + gr, uk_map + gpr, -chi_p*nup_sigf[gp] * volume);
          }

          //===== Delayed
          if (not lag_precursors)
          {
            const double* chi_d = xs->chi_delayed[g].data();
            const double* nud_sigf = xs->nu_delayed_sigma_f.data();
            const double* lambda = xs->precursor_lambda.data();
            const double* gamma = xs->precursor_yield.data();

            /* This section of code computes the matrix coefficient
             * corresponding to the precursor substitution term arising when
             * precursors are treated implicitly. This term is similar to the
             * normal fission term with the exception of having an inner sum
             * over all precursor species. The first loop computes sum over
             * precursor species. This provides a factor which is multiplied
             * by nu_delayed_sigma_f when contributing to the matrix. The
             * coefficient is computed ahead of time for the sake of efficiency.
             * Using nested loops would require either this factor to be
             * computed \p n_gsg times or for many more expensive sparse matrix
             * add operations. */
            double coeff= 0.0;
            for (size_t j = 0; j < xs->n_precursors; ++j)
              coeff +=
                chi_d[j]*lambda[j]/(1.0 + eff_dt*lambda[j]) * gamma[j]*eff_dt;

            // Contribute the delayed fission term to the matrix.
            for (size_t gpr = 0; gpr < n_gsg; ++gpr)
            {
              const size_t gp = groupset.groups[gpr];
              A.add(uk_map + gr, uk_map + gpr, coeff*nud_sigf[gp] * volume);
            }
          }//if not lag precursors
        }//if prompt + delayed fission
      }//if fissile
    }//for group

    //============================== Loop over faces
    for (const auto& face : cell.faces)
    {
      //==================== Interior faces
      if (face.has_neighbor)
      {
        // Get neighbor info
        const auto& nbr_cell = mesh->cells[face.neighbor_id];
        const auto nbr_xs_id = matid_to_xs_map[nbr_cell.material_id];
        const auto& nbr_xs = material_xs[nbr_xs_id];

        const size_t nbr_uk_map = nbr_cell.id * n_gsg;

        // Geometric data
        const double d_pf = cell.centroid.distance(face.centroid);
        const double d_pn = cell.centroid.distance(nbr_cell.centroid);
        const double w = d_pf / d_pn; // harmonic mean weight factor

        // Diffusion coefficients
        const double* D = xs->diffusion_coeff.data();
        const double* D_nbr = nbr_xs->diffusion_coeff.data();

        for (size_t gr = 0; gr < n_gsg; ++gr)
        {
          const size_t g = groupset.groups[gr];
          const auto D_eff = 1.0/(w/D[g] + (1.0 - w)/D_nbr[g]);

          A.add(uk_map + gr, uk_map + gr, D_eff/d_pn * face.area);
          A.add(uk_map + gr, nbr_uk_map + gr, -D_eff/d_pn * face.area);
        }
      }//if interior face

      //==================== Boundary faces
      else
      {
        const auto bndry_id = face.neighbor_id;
        const auto bndry_type = boundary_info[bndry_id].first;

        //========== Dirichlet boundary term
        if (bndry_type == BoundaryType::ZERO_FLUX ||
            bndry_type == BoundaryType::DIRICHLET)
        {
          const double* D = xs->diffusion_coeff.data();
          const double d_pf = cell.centroid.distance(face.centroid);

          for (size_t gr = 0; gr < n_gsg; ++gr)
          {
            const size_t g = groupset.groups[gr];
            A.add(uk_map + gr, uk_map + gr, D[g]/d_pf * face.area);
          }
        }//if Dirichlet

        //========== Robin boundary term
        else if (bndry_type == BoundaryType::VACUUM ||
                 bndry_type == BoundaryType::MARSHAK ||
                 bndry_type == BoundaryType::ROBIN)
        {
          const double* D = xs->diffusion_coeff.data();
          const double d_pf = cell.centroid.distance(face.centroid);

          for (size_t gr = 0; gr < n_gsg; ++gr)
          {
            const size_t g = groupset.groups[gr];

            const auto& bndry = boundaries[bndry_id][g];
            const auto bc = std::static_pointer_cast<RobinBoundary>(bndry);

            double val = bc->a*D[g]/(bc->b*D[g] + bc->a*d_pf);
            A.add(uk_map + gr, uk_map + gr, val * face.area);
          }
        }//if Robin
      }//if boundary face
    }//for face
  }//for cell
}


void
TransientSolver::assemble_matrices()
{
  for (auto& groupset : groupsets)
  {
    if (solution_technique == SolutionTechnique::GROUPSET_WISE)
      assemble_transient_matrix(groupset, NO_ASSEMBLER_FLAGS);
    else
      assemble_transient_matrix(groupset,
                                ASSEMBLE_SCATTER | ASSEMBLE_FISSION);

    linear_solver->set_matrix(groupset.A);
  }
}
