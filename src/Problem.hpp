#ifndef NAVIER_HPP
#define NAVIER_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_in.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include "PreconditionBlockTriangular.hpp"
#include "SimulationSettings.hpp"

#include <fstream>
#include <iostream>

using namespace dealii;

// Class implementing a solver for the NavierStokes-Stokes problem.
class NavierStokes
{
public:
    // Physical dimension (1D, 2D, 3D)
    static constexpr unsigned int dim = 3;

    // Function for the forcing term.
    class ForcingTerm : public Function<dim>
    {
        // get_time, set_time, advance_time
    public:
        virtual void
        vector_value(const Point<dim> & /*p*/,
                     Vector<double> &values) const override
        {
            for (unsigned int i = 0; i < dim - 1; ++i)
                values[i] = 0.0;

            values[dim - 1] = -g;
        }

        virtual double
        value(const Point<dim> & /*p*/,
              const unsigned int component = 0) const override
        {
            if (component == dim - 1)
                return -g;
            else
                return 0.0;
        }

    protected:
        const double g = 0.0;
    };

    // Function for inlet velocity. This actually returns an object with four
    // components (one for each velocity component, and one for the pressure), but
    // then only the first three are really used (see the component mask when
    // applying boundary conditions at the end of assembly). If we only return
    // three components, however, we may get an error message due to this function
    // being incompatible with the finite element space.
    class InletVelocity : public Function<dim>
    {
    public:
        InletVelocity(NavierStokes &parent)
            : Function<dim>(dim + 1), parent(parent)
        {
        }

        virtual void
        vector_value(const Point<dim> /*&p*/, Vector<double> &values) const
        {
            // values[0] = -alpha * p[1] * (2.0 - p[1]) * (1.0 - p[2]) * (2.0 - p[2]);
            values[3] = 0.0; // Unused component
            for (unsigned int i = 0; i < dim; ++i)
                values[i] = parent.get_current_inlet_velocity(i);
        }

        virtual double
        value(const Point<dim> /*&p*/, const unsigned int component = 0) const
        {
            if (component == 3) // Unused component
                // return -alpha * p[1] * (2.0 - p[1]) * (1.0 - p[2]) * (2.0 - p[2]);
                return 0.0;
            else
            {
                return parent.get_current_inlet_velocity(component);
            }
        }

    protected:
        NavierStokes &parent;
    };

    // Since we're working with block matrices, we need to make our own
    // preconditioner class. A preconditioner class can be any class that exposes
    // a vmult method that applies the inverse of the preconditioner.

    // Identity preconditioner.
    class PreconditionIdentity
    {
    public:
        // Application of the preconditioner: we just copy the input vector (src)
        // into the output vector (dst).
        void
        vmult(TrilinosWrappers::MPI::BlockVector &dst,
              const TrilinosWrappers::MPI::BlockVector &src) const
        {
            dst = src;
        }

    protected:
    };

    // Constructor.
    NavierStokes(ns_sim_settings::SimulationSettings &settings)
        : settings(settings), mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)), mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)), pcout(std::cout, mpi_rank == 0), inlet_velocity(*this), mesh(MPI_COMM_WORLD)
    {
    }

    // Setup system.
    void
    setup();

    // Assemble system for each iteration of the Newton method. We also assemble the pressure mass matrix (needed for the
    // preconditioner).
    void
    assemble_system(const bool initial_step);

    // Solve Newton step.
    void
    solve_newton_step();

    // Solve Newton method.
    void
    solve_newton();

    // Output results.
    void
    output() const;

    void increment_time_step();

    double get_current_coeff_nu() const;

    double get_current_inlet_velocity(const unsigned int component = 0) const;

protected:
    // MPI parallel. /////////////////////////////////////////////////////////////

    unsigned long current_time_step = 0;

    ns_sim_settings::SimulationSettings& settings;

    // Number of MPI processes.
    const unsigned int mpi_size;

    // This MPI process.
    const unsigned int mpi_rank;

    // Parallel output stream.
    ConditionalOStream pcout;

    // Problem definition. ///////////////////////////////////////////////////////

    // Kinematic viscosity [m2/s].
    // const double nu = 1;

    // Gamma parameter
    // const double gamma = 1.0;

    // Outlet pressure [Pa].
    // const double p_out = 10;

    // Forcing term.
    ForcingTerm forcing_term;

    // Inlet velocity.
    InletVelocity inlet_velocity;

    // Discretization. ///////////////////////////////////////////////////////////

    // Mesh refinement.
    // const unsigned int N;

    // Polynomial degree used for velocity.
    // const unsigned int degree_velocity;

    // Polynomial degree used for pressure.
    // const unsigned int degree_pressure;

    // Mesh.
    parallel::fullydistributed::Triangulation<dim> mesh;

    // Finite element space.
    std::unique_ptr<FiniteElement<dim>> fe;

    // Quadrature formula.
    std::unique_ptr<Quadrature<dim>> quadrature;

    // Quadrature formula for face integrals.
    std::unique_ptr<Quadrature<dim - 1>> quadrature_face;

    // DoF handler.
    DoFHandler<dim> dof_handler;

    // DoFs owned by current process.
    IndexSet locally_owned_dofs;

    // DoFs owned by current process in the velocity and pressure blocks.
    std::vector<IndexSet> block_owned_dofs;

    // DoFs relevant to the current process (including ghost DoFs).
    IndexSet locally_relevant_dofs;

    // DoFs relevant to current process in the velocity and pressure blocks.
    std::vector<IndexSet> block_relevant_dofs;

    // System matrix.
    TrilinosWrappers::BlockSparseMatrix jacobian_matrix;

    // Pressure mass matrix, needed for preconditioning. We use a block matrix for
    // convenience, but in practice we only look at the pressure-pressure block.
    TrilinosWrappers::BlockSparseMatrix pressure_mass;

    // Right-hand side vector in the linear system.
    TrilinosWrappers::MPI::BlockVector residual_vector;

    // Solutin increment (without ghost elements).
    TrilinosWrappers::MPI::BlockVector delta_owned;

    // System solution (without ghost elements).
    TrilinosWrappers::MPI::BlockVector solution_owned;

    // System solution (including ghost elements).
    TrilinosWrappers::MPI::BlockVector solution;
};

#endif