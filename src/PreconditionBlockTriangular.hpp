#ifndef PRECONDITION_BLOCK_TRIANGULAR_HPP
#define PRECONDITION_BLOCK_TRIANGULAR_HPP

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

using namespace dealii;

// Implementation of the preconditioner as described in the
// deal.ii step 57 example: the approach used is taken from:
// "An Augmented Lagrangian-Based Approach to the Oseen Problem", M. Benzi and M. Olshanskii, SIAM J. SCI. COMPUT. 2006.
//
// An easy-to-read explanation at this web page: https://epubs.siam.org/doi/epdf/10.1137/18M1219370
// The idea is basically always the same: using a Block Proconditioner to solve the saddle point problem
// as an attempt to generate solvers whose convergence is indipendent of parametric values
// and mesh size h.
// In this case, the use of Block Preconditioners allows to deal with the problem using the block gaussian elimitation:
// this implies reducing the problem considered from solving the coupled linear system to solving smaller separate linear systems
// involving the A matrix and the S = - B * A^-1 * B^T (Shur complement) matrix.
// As a note: this approach is totally anologous to the one followed by the lab-09. Yet in the Stokes equation case
// the Shur complement was simply spectrally equivalent to the Pressure Mass Matrix weighted by the viscosity coefficient.
// Still this method doesn't really fare well in the case the Reynolds number gets too big.
// The method proposed by Benzi and Olshanskii avoids this problem by adding a term in the equation that does not change the
// continous solution, while at the same time allows for modification of the Shur complement.

// Consideriamo Q * P^-1 * y = b con Q la matrice ottenuta da saddle point prob, con l'agggiunta del termine
// gamma * B^T * W^-1 * B nel blocco in alto a sinistra: in questo caso si tratta del termine aggiuntivo
// coerente alla selta del precondizionatore
/* NAVIER - STOKES
class PreconditionBlockTriangular
{
public:
  // Initialize the preconditioner, given the velocity stiffness matrix, the
  // pressure mass matrix.
  void
  initialize(const TrilinosWrappers::SparseMatrix &velocity_stiffness_,
             const TrilinosWrappers::SparseMatrix &pressure_mass_,
             const TrilinosWrappers::SparseMatrix &B_T_,
             const double &nu_,
             const double & gamma_) // jacobian_matrix.block(0, 1)
  {
    velocity_stiffness = &velocity_stiffness_;
    pressure_mass = &pressure_mass_;
    B_T = &B_T_;
    nu = nu_;
    gamma = gamma_;

    preconditioner_velocity.initialize(velocity_stiffness_);
    preconditioner_pressure.initialize(pressure_mass_);
  }

  // Application of the preconditioner.
  // P * x = y  --> G * P^-1 * y = b con P^-1 * y = x
  void
  vmult(TrilinosWrappers::MPI::BlockVector &dst,
        const TrilinosWrappers::MPI::BlockVector &src) const
  {

    {

      // Application of the lower half of the preconditioner
      SolverControl solver_control_pressure(1000,
                                            1e-2 * src.block(1).l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_pressure(
          solver_control_pressure);
      solver_cg_pressure.solve(*pressure_mass,
                               dst.block(1),
                               src.block(1),
                               preconditioner_pressure);
      dst.block(1) *= -(nu + gamma);
    }

    {

      // Application of the upper half of the preconditioner
      // "using the data for the lower half"
      tmp.reinit(src.block(0));
      B_T->vmult(tmp, dst.block(1));
      tmp.sadd(-1.0, src.block(0));

      // SolverControl solver_control_velocity(1000,
      //                                       1e-2 * src.block(0).l2_norm());
      // SolverGMRES<TrilinosWrappers::MPI::Vector> solver_cg_velocity(
      //     solver_control_velocity);
      // solver_cg_velocity.solve(*velocity_stiffness,
      //                           dst.block(0),
      //                           tmp,
      //                           preconditioner_velocity);

      velocity_stiffness->vmult(dst.block(0), tmp);

      
    }
  }

  

    

protected:
  // Viscosity coefficient
  double nu;
  // Gamma coefficient
  double gamma;
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
*/
 // SOLO STOKES
class PreconditionBlockTriangular
  {
  public:
    // Initialize the preconditioner, given the velocity stiffness matrix, the
    // pressure mass matrix.
    void
    initialize(const TrilinosWrappers::SparseMatrix &velocity_stiffness_,
               const TrilinosWrappers::SparseMatrix &pressure_mass_,
               const TrilinosWrappers::SparseMatrix &B_,
               const double &nu_,
             const double & gamma_)
    {
      velocity_stiffness = &velocity_stiffness_;
      pressure_mass = &pressure_mass_;
      B = &B_;

      preconditioner_velocity.initialize(velocity_stiffness_);
      preconditioner_pressure.initialize(pressure_mass_);
    }

    // Application of the preconditioner.
    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      SolverControl solver_control_velocity(1000,
                                            1e-2 * src.block(0).l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_velocity(
          solver_control_velocity);
      solver_cg_velocity.solve(*velocity_stiffness,
                               dst.block(0),
                               src.block(0),
                               preconditioner_velocity);

      tmp.reinit(src.block(1));
      B->vmult(tmp, dst.block(0));
      tmp.sadd(-1.0, src.block(1));

      SolverControl solver_control_pressure(1000,
                                            1e-2 * src.block(1).l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_pressure(
          solver_control_pressure);
      solver_cg_pressure.solve(*pressure_mass,
                               dst.block(1),
                               tmp,
                               preconditioner_pressure);
    }

  protected:
    // Velocity stiffness matrix.
    const TrilinosWrappers::SparseMatrix *velocity_stiffness;

    // Preconditioner used for the velocity block.
    TrilinosWrappers::PreconditionILU preconditioner_velocity;

    // Pressure mass matrix.
    const TrilinosWrappers::SparseMatrix *pressure_mass;

    // Preconditioner used for the pressure block.
    TrilinosWrappers::PreconditionILU preconditioner_pressure;

    // B matrix.
    const TrilinosWrappers::SparseMatrix *B;

    // Temporary vector.
    mutable TrilinosWrappers::MPI::Vector tmp;
  };

#endif