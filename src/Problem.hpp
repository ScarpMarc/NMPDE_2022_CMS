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

// #include "PreconditionBlockTriangular.hpp"
#include "SimulationSettings.hpp"
#include "PreconditionSIMPLE.hpp"

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

    // The initial solution to set up the system.
    class InitialSolution : public Function<dim>
    {
    public:
        InitialSolution()
            : Function<dim>(dim + 1)
        {
        }

        virtual void
        vector_value(const Point<dim> & /*p*/, Vector<double> &values) const override
        {
            for (unsigned int i = 0; i < dim + 1; ++i)
                values[i] = 0.0;
        }

        virtual double
        value(const Point<dim> & /*p*/, const unsigned int /*component*/ = 0) const override
        {
            return 0.0;
        }
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
        InletVelocity(const std::array<double, 3> &inlet_velocity_start, const std::array<double, 3> &inlet_velocity_end, const unsigned int &time_steps_pre_ramp, const unsigned int &time_steps_ramp, const unsigned int &time_steps_per_second)
            : Function<dim>(dim + 1), velocity_start(inlet_velocity_start.begin(), inlet_velocity_start.end()), 
            velocity_end(inlet_velocity_end.begin(), inlet_velocity_end.end()), time_steps_pre_ramp(time_steps_pre_ramp), time_steps_ramp(time_steps_ramp), time_steps_per_second(time_steps_per_second)
        {
        }

        /*
        if (current_time_step < settings.get_time_steps_pre_ramp())
    return settings.get_inlet_velocity_start()[component];
  else if (current_time_step >= settings.get_time_steps_pre_ramp() + settings.get_time_steps_ramp())
    return settings.get_inlet_velocity_end()[component];
  else
    return settings.get_inlet_velocity_start()[component] +
           (((double)current_time_step) / (double)settings.get_time_steps_per_second()) *
               (settings.get_inlet_velocity_end()[component] - settings.get_inlet_velocity_start()[component]);


        */

        virtual void
        vector_value(const Point<dim> & /*p*/, Vector<double> &values) const override
        {
            // values[0] = -alpha * p[1] * (2.0 - p[1]) * (1.0 - p[2]) * (2.0 - p[2]);
            values[3] = 0.0; // Unused component
            for(unsigned int i = 0; i < dim; ++i)
                values[i] = velocity_start[i] + (velocity_end[i] - velocity_start[i]) * get_time();
        }

        virtual double
        value(const Point<dim> & /*p*/, const unsigned int component = 0) const override
        {
            if (component == 3) // Unused component
                // return -alpha * p[1] * (2.0 - p[1]) * (1.0 - p[2]) * (2.0 - p[2]);
                return 0.0;
            else
            {
                return velocity_start[component] + (velocity_end[component] - velocity_start[component]) * get_time();
            }
        }

        /**
         * @brief Set the time step.
         * @details The velocity is kept at the initial value before the ramp, so we set the time to 0. Then it increases linearly during the ramp, and is set to the final value after the ramp. This is done so that the value and vector_value functions need to just return (v_end - v_start) at time t, without any if statements.
         */
        void set_time_step(const unsigned long &time_step)
        {
            if (time_step < time_steps_pre_ramp)
                set_time(0.0);
            else if (time_step >= time_steps_pre_ramp + time_steps_ramp)
                set_time((double)(time_steps_pre_ramp + time_steps_ramp) / (double)time_steps_per_second);
            else
                set_time((double)time_step / (double)time_steps_per_second);
        }

    protected:
        // const NavierStokes &parent;

        /*inline constexpr double get_cur_velocity_value(const unsigned int component = 0)
        {
            if (get_time() < ((double)time_steps_pre_ramp / (double)time_steps_per_second))
                return inlet_velocity_start[component];
            else if (get_time() >= (((double)time_steps_pre_ramp + (double)time_steps_ramp) / (double)time_steps_per_second))
                return inlet_velocity_end[component];
            else
                return inlet_velocity_start[component] +
                       get_time() *
                           (inlet_velocity_end[component] - inlet_velocity_start[component]);
        }*/

        const Vector<double> velocity_start;
        const Vector<double> velocity_end;
        const unsigned int time_steps_pre_ramp;
        const unsigned int time_steps_ramp;
        const unsigned int time_steps_per_second;
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

    /*
    const std::vector<unsigned int> surfaces_walls;

    const std::vector<unsigned int> surfaces_inlets;

    const std::vector<unsigned int> surfaces_outlets;

    const std::vector<unsigned int> rn surfaces_free_slip;

    */

    // Constructor.
    NavierStokes(ns_sim_settings::SimulationSettings &settings)
        : settings(settings), mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)), mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)), pcout(std::cout, mpi_rank == 0),
          surfaces_walls(settings.get_surfaces_walls()), surfaces_inlets(settings.get_surfaces_inlets()), surfaces_outlets(settings.get_surfaces_outlets()), surfaces_free_slip(settings.get_surfaces_free_slip()),
          inlet_velocity(settings.get_inlet_velocity_start(), settings.get_inlet_velocity_end(), settings.get_time_steps_pre_ramp(), settings.get_time_steps_ramp(), settings.get_time_steps_per_second()), mesh(MPI_COMM_WORLD)
    {
    }

    // Setup system.
    void setup();

    /**
     * @brief Assemble the system matrix and right-hand side.
     */
    void
        assemble_system(/*const AssemblyType &type*/);

    // Solve Newton step.
    // void
    // solve_newton_step();

    // Solve Newton method.
    // void
    // solve_newton();

    /**
     * @brief The entry point of the solving process.
     * @details Evaluates Reynold's number and decides whether to employ the continuation method at the start or not. Then, it calls the solve_time_step() function for each time step, and increments the time step counter after each iteration.
     */
    void solve();

    /**
     * @brief Solves the system for the current time step.
     * @details This function also exports the solution to a file. Note that this function does not increment the time step counter, nor does it perform any checks on Reynold's number.
     */
    void solve_time_step();

    /**
     * @brief Computes the initial guess by using the continuation method.
     * @details The function will iterate through the values of the kinematic viscosity and solve the system for each value. The solution of the previous iteration will be used as the initial guess for the next iteration. Then, the solution of the last iteration will be used as the initial guess for the actual simulation.
     * @param nu_start Starting value of the kinematic viscosity.
     * @param nu_end Ending value of the kinematic viscosity.
     * @param num_steps Number of steps to take between the starting and ending values.
     */
    void compute_initial_guess(const double &nu_start, const double &nu_end, const unsigned int &num_steps);

    /**
     * @brief Estimate Reynold's number.
     * @return double The estimated Reynold's number.
     */
    double estimate_reynolds_number() const;

    // Output results.
    /**
     * @brief Outputs the solution to a file, at the current time step.
     * @details Called at the end of each time step.
     */
    void output() const;

    void increment_time_step();

    double get_current_coeff_nu() const;

    /*
        double get_current_inlet_velocity(const unsigned int &component = 0) const;
        void get_current_inlet_velocity(Vector<double> &output) const;
    */

    // enum AssemblyType
    //{
    //     STOKES_ONLY,
    //     NAVIER_STOKES
    // };

protected:
    // MPI parallel. /////////////////////////////////////////////////////////////

    ns_sim_settings::SimulationSettings &settings;

    unsigned long current_time_step = 0;

    // Number of MPI processes.
    const unsigned int mpi_size;

    // This MPI process.
    const unsigned int mpi_rank;

    // Parallel output stream.
    ConditionalOStream pcout;

    // Problem definition. ///////////////////////////////////////////////////////

    const std::vector<unsigned int> surfaces_walls;

    const std::vector<unsigned int> surfaces_inlets;

    const std::vector<unsigned int> surfaces_outlets;

    const std::vector<unsigned int> surfaces_free_slip;

    // Kinematic viscosity [m2/s].
    // const double nu = 1;

    // Gamma parameter
    // const double gamma = 1.0;

    // Outlet pressure [Pa].
    // const double p_out = 10;

    const double max_reynolds_number_before_continuation = 500;

    // Forcing term.
    ForcingTerm forcing_term;

    // Inlet velocity.
    InletVelocity inlet_velocity;

    InitialSolution initial_solution;

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
    // TrilinosWrappers::BlockSparseMatrix pressure_mass;

    // Diagonal Inverse of F matrix, needed for preconditioning. We use a blockvector for convenience.
    TrilinosWrappers::MPI::BlockVector D_inv;

    // Right-hand side vector in the linear system.
    TrilinosWrappers::MPI::BlockVector residual_vector;

    // Solution increment (without ghost elements).
    // TrilinosWrappers::MPI::BlockVector delta_owned;

    // System solution (without ghost elements).
    TrilinosWrappers::MPI::BlockVector solution_owned;

    // System solution (including ghost elements).
    TrilinosWrappers::MPI::BlockVector solution;

    // System solution at previous time step.
    TrilinosWrappers::MPI::BlockVector solution_old;
};

#endif