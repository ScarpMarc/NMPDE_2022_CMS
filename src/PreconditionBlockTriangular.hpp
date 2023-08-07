// Block-triangular preconditioner. (Block Shur Preconditioner, as described in https://www.dealii.org/current/doxygen/deal.II/step_57.html)
class PreconditionBlockTriangular
{
public:
  // Initialize the preconditioner, given:
  // - the velocity stiffness matrix
  // - the pressure mass matrix
  // - the B matrix
  // - the viscosity
  // - the gamma parameter
  void
  initialize(const double &gamma_,
             const double &viscosity_,
             const TrilinosWrappers::SparseMatrix &velocity_stiffness_modified_,
             const TrilinosWrappers::SparseMatrix &pressure_mass_,
             const TrilinosWrappers::SparseMatrix &B_T_)
  {
    gamma = &gamma_;
    viscosity = &viscosity_;
    velocity_stiffness = &velocity_stiffness_modified_;
    pressure_mass = &pressure_mass_;
    B_T = &B_T_;

    preconditioner_velocity.initialize(velocity_stiffness_modified_);
    preconditioner_pressure.initialize(pressure_mass_);
  }

  // Application of the preconditioner.
  void
  vmult(TrilinosWrappers::MPI::BlockVector &dst,
        const TrilinosWrappers::MPI::BlockVector &src) const
  {

    // Effect lower part of the preconditioner (Shur complement)
    // Computing dst_1 = -(gamma + viscosity)*Mp^-1 * src_1
    {
      SolverControl solver_control_pressure(1000,
                                            1e-2 * src.block(1).l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_pressure(
          solver_control_pressure);
      solver_cg_pressure.solve(*pressure_mass,
                               dst.block(1),
                               src.block(1),
                               preconditioner_pressure);
      dst.block(1) *= -(*gamma + *viscosity);
    }

    // Effect of the upper part of the block-triangular preconditioner
    // Computing dest_1 = A^-1 * (src_0 - B^T * dst_1)

    {
      tmp.reinit(src.block(0));
      B_T->vmult(tmp, dst.block(1));
      tmp.sadd(-1.0, src.block(0));

      // Change the iterative solver for the velocity sub-block (A is not simmetric anymore)
      SolverControl solver_control_velocity(1000,
                                            1e-2 * src.block(0).l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres_velocity(
          solver_control_velocity);
      solver_gmres_velocity.solve(*velocity_stiffness,
                                  dst.block(0),
                                  tmp,
                                  preconditioner_velocity);
    }
  }

protected:
  // Gamma parameter.
  const double *gamma;

  // Viscosity.
  const double *viscosity;

  // Velocity stiffness matrix.
  const TrilinosWrappers::SparseMatrix *velocity_stiffness;

  // Preconditioner used for the velocity block.
  TrilinosWrappers::PreconditionILU preconditioner_velocity;

  // Pressure mass matrix.
  const TrilinosWrappers::SparseMatrix *pressure_mass;

  // Preconditioner used for the pressure block.
  TrilinosWrappers::PreconditionILU preconditioner_pressure;

  // B matrix.
  const TrilinosWrappers::SparseMatrix *B_T;

  // Temporary vector.
  mutable TrilinosWrappers::MPI::Vector tmp;
};