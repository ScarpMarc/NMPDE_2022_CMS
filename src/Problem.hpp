

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
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <array>
#include <iostream>

#include "SimulationSettings.hpp"

using namespace dealii;

class Problem
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

  // Function for inlet velocity. This actually returns an object with four
  // components (one for each velocity component, and one for the pressure), but
  // then only the first three are really used (see the component mask when
  // applying boundary conditions at the end of assembly). If we only return
  // three components, however, we may get an error message due to this function
  // being incompatible with the finite element space.
  class InletVelocity : public Function<dim>
  {
  public:
    InletVelocity(const ns_sim_settings::SimulationSettings &simulation_settings_)
        : Function<dim>(dim + 1), components(simulation_settings_.inlet_velocity)
    {
    }

    virtual void
    vector_value(const Point<dim> & /*p*/, Vector<double> &values) const override
    {

      for (unsigned int i = 0; i < dim; ++i)
        values[i] = components[i] * factor;

      values[3] = 0.0;
    }

    virtual double
    value(const Point<dim> & /*p*/, const unsigned int component = 0) const override
    {
      if (component == 3)
        return 0.0;
      else
        return components[component] * factor;
    }

    void set_factor(const double &factor) { factor = factor; }
    double get_factor() const { return factor; }

  protected:
    const std::array<double, 3> components;

    double factor = 0;
  };

  // Constructor.
  Problem(const ns_sim_settings::SimulationSettings &simulation_settings)
      : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
        mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
        simulation_settings(simulation_settings),
        pcout(std::cout, mpi_rank == 0),
        inlet_velocity(simulation_settings),
        /*file_name(SimulationSettings.file_name), degree_velocity(SimulationSettings.degree_velocity), degree_pressure(SimulationSettings.degree_pressure), nu(SimulationSettings.coeff_nu), p_out(SimulationSettings.outlet_pressure), */ mesh(MPI_COMM_WORLD)
  {
  }

  // Setup system.
  void setup();

  void solve_time_step();

  // Assemble system. We also assemble the pressure mass matrix (needed for the
  // preconditioner).
  void assemble(bool only_rhs = false, double nonlinearita = 1.0);

  void Problem::compute_initial_guess(double step_size);

  // Solve linear system.
  void solveLinearSystem(); // risolve il sistema lineare

  // solve Newton Method
  void solveNewtonMethod(); // risolve il problema del metodo di newthon
                            // chiamando piu volte solveLinearSystem
  // solve problem in time
  // void                          // risolve il problema per ogni istante di tempo
  // solve();                      // TODO

  // Output results.
  void output(const unsigned long &time_step);

protected:
  // MPI parallel. /////////////////////////////////////////////////////////////
  unsigned long current_time;

  ns_sim_settings::SimulationSettings simulation_settings;

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Problem definition. ///////////////////////////////////////////////////////

  // Kinematic viscosity [m2/s].
  // const double nu = 1.;

  // Gamma parameter
  // const double gamma = 1.;

  // Outlet pressure [Pa].
  // const double p_out = 10.;

  // Forcing term.
  ForcingTerm forcing_term;

  // Inlet velocity.
  InletVelocity inlet_velocity;

  // Discretization. ///////////////////////////////////////////////////////////

  // Mesh refinement.
  // std::string filename;

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

  // Pressure mass matrix, needed for preconditioning. We use a block matrix for
  // convenience, but in practice we only look at the pressure-pressure block.
  TrilinosWrappers::BlockSparseMatrix pressure_mass;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::BlockVector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::BlockVector solution;

  //
  //  variabili relative alla linearizzazione del sistema =>
  //

  // Jacobian matrix.
  TrilinosWrappers::BlockSparseMatrix jacobian_matrix;

  // Residual vector.
  TrilinosWrappers::MPI::BlockVector residual_vector;

  // Increment of the solution between Newton iterations.
  TrilinosWrappers::MPI::BlockVector delta_owned;
};
