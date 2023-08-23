#include "Problem.hpp"
#include <deal.II/base/timer.h>
#include <algorithm>

double NavierStokes::get_current_coeff_nu() const
{
  // For now
  return settings.get_coeff_nu();
}

/*double NavierStokes::get_current_inlet_velocity(const unsigned int &component) const
{

}

void NavierStokes::get_current_inlet_velocity(Vector<double> &output) const
{
  for (unsigned int i = 0; i < dim; ++i)
  {
    output[i] = get_current_inlet_velocity(i);
  }
}*/

void NavierStokes::increment_time_step()
{
  ++current_time_step;
  inlet_velocity.set_time_step(current_time_step);
  forcing_term.set_time((double)current_time_step / (double)settings.get_time_steps_per_second());
  initial_solution.set_time((double)current_time_step / (double)settings.get_time_steps_per_second());
}

inline double NavierStokes::estimate_reynolds_number() const
{
  Point<dim> temp;
  Vector<double> cur_velocity(dim + 1);
  inlet_velocity.vector_value(temp, cur_velocity);
  return (settings.get_coeff_rho() * cur_velocity.l2_norm() * settings.get_characteristic_length()) / get_current_coeff_nu();
}

void NavierStokes::setup()
{
  // Create the mesh.
  {
    pcout << "Initializing the mesh" << std::endl;

    Triangulation<dim> mesh_serial;

    GridIn<dim> grid_in;
    grid_in.attach_triangulation(mesh_serial);

    const std::string mesh_file_name =
        settings.get_file_name();

    std::ifstream grid_in_file(mesh_file_name);
    grid_in.read_msh(grid_in_file);

    GridTools::partition_triangulation(mpi_size, mesh_serial);
    const auto construction_data = TriangulationDescription::Utilities::
        create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
    mesh.create_triangulation(construction_data);

    pcout << "  Number of elements = " << mesh.n_global_active_cells()
          << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);

    for (unsigned int i = 0; i < mpi_size; ++i)
    {
      if (i == mpi_rank)
      {
        cout << "    Number of elements for processor " << i << " = " << mesh.n_active_cells()
             << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space.
  {
    pcout << "Initializing the finite element space" << std::endl;

    const FE_SimplexP<dim> fe_scalar_velocity(settings.get_degree_velocity());
    const FE_SimplexP<dim> fe_scalar_pressure(settings.get_degree_pressure());
    fe = std::make_unique<FESystem<dim>>(fe_scalar_velocity,
                                         dim,
                                         fe_scalar_pressure,
                                         1);

    pcout << "  Velocity degree:           = " << fe_scalar_velocity.degree
          << std::endl;
    pcout << "  Pressure degree:           = " << fe_scalar_pressure.degree
          << std::endl;
    pcout << "  DoFs per cell              = " << fe->dofs_per_cell
          << std::endl;

    quadrature = std::make_unique<QGaussSimplex<dim>>(fe->degree + 1);

    pcout << "  Quadrature points per cell = " << quadrature->size()
          << std::endl;

    quadrature_face = std::make_unique<QGaussSimplex<dim - 1>>(fe->degree + 1);

    pcout << "  Quadrature points per face = " << quadrature_face->size()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  {
    pcout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    // We want to reorder DoFs so that all velocity DoFs come first, and then
    // all pressure DoFs.
    std::vector<unsigned int> block_component(dim + 1, 0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, block_component);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    // Besides the locally owned and locally relevant indices for the whole
    // system (velocity and pressure), we will also need those for the
    // individual velocity and pressure blocks.
    std::vector<types::global_dof_index> dofs_per_block =
        DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
    const unsigned int n_u = dofs_per_block[0];
    const unsigned int n_p = dofs_per_block[1];

    block_owned_dofs.resize(2);
    block_relevant_dofs.resize(2);
    block_owned_dofs[0] = locally_owned_dofs.get_view(0, n_u);
    block_owned_dofs[1] = locally_owned_dofs.get_view(n_u, n_u + n_p);
    block_relevant_dofs[0] = locally_relevant_dofs.get_view(0, n_u);
    block_relevant_dofs[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p);

    pcout << "  Number of DoFs: " << std::endl;
    pcout << "    velocity = " << n_u << std::endl;
    pcout << "    pressure = " << n_p << std::endl;
    pcout << "    total    = " << n_u + n_p << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    pcout << "Initializing the linear system" << std::endl;

    pcout << "  Initializing the sparsity pattern" << std::endl;

    // Velocity DoFs interact with other velocity DoFs (the weak formulation has
    // terms involving u times v), and pressure DoFs interact with velocity DoFs
    // (there are terms involving p times v or u times q). However, pressure
    // DoFs do not interact with other pressure DoFs (there are no terms
    // involving p times q). We build a table to store this information, so that
    // the sparsity pattern can be built accordingly.
    Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
    for (unsigned int c = 0; c < dim + 1; ++c)
    {
      for (unsigned int d = 0; d < dim + 1; ++d)
      {
        if (c == dim && d == dim) // pressure-pressure term
          coupling[c][d] = DoFTools::none;
        else // other combinations
          coupling[c][d] = DoFTools::always;
      }
    }

    TrilinosWrappers::BlockSparsityPattern sparsity(block_owned_dofs,
                                                    MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, coupling, sparsity);
    sparsity.compress();

    // We also build a sparsity pattern for the pressure mass matrix.
    /*for (unsigned int c = 0; c < dim + 1; ++c)
    {
      for (unsigned int d = 0; d < dim + 1; ++d)
      {
        if (c == dim && d == dim) // pressure-pressure term
          coupling[c][d] = DoFTools::always;
        else // other combinations
          coupling[c][d] = DoFTools::none;
      }
    }
    TrilinosWrappers::BlockSparsityPattern sparsity_pressure_mass(
        block_owned_dofs, MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler,
                                    coupling,
                                    sparsity_pressure_mass);
    sparsity_pressure_mass.compress();*/

    pcout << "  Initializing the matrix" << std::endl;
    jacobian_matrix.reinit(sparsity);
    // pressure_mass.reinit(sparsity_pressure_mass);

    pcout << "  Initializing the system right-hand side" << std::endl;
    residual_vector.reinit(block_owned_dofs, MPI_COMM_WORLD);
    D_inv.reinit(block_owned_dofs, MPI_COMM_WORLD);
    pcout << "  Initializing the solution vector" << std::endl;
    solution_owned.reinit(block_owned_dofs, MPI_COMM_WORLD);
    // delta_owned.reinit(block_owned_dofs, MPI_COMM_WORLD);
    solution.reinit(block_owned_dofs, block_relevant_dofs, MPI_COMM_WORLD);
    solution_old = solution;
  }
}

void NavierStokes::assemble_system(/*const AssemblyType &type*/)
{
  // Initial step: stokes only. Needed
  pcout << "===============================================" << std::endl;
  pcout << "Assembling the system";
  // if (type == STOKES_ONLY)
  //   pcout << " (Stokes only)";
  pcout << std::endl;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q = quadrature->size();
  const unsigned int n_q_face = quadrature_face->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
  FEFaceValues<dim> fe_face_values(*fe,
                                   *quadrature_face,
                                   update_values | update_normal_vectors |
                                       update_JxW_values);

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);
  Vector<double> cell_diag(dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  jacobian_matrix = 0.0;
  residual_vector = 0.0;
  D_inv = 0.0;
  // pressure_mass = 0.0;

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  // We use this vector to store the solution at the previous time step.
  std::vector<Tensor<1, dim>> velocity_old_loc(n_q);

  // Inlet velocity and forcing term times are already set.

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (!cell->is_locally_owned())
      continue;

    fe_values.reinit(cell);

    cell_matrix = 0.0;
    cell_rhs = 0.0;
    cell_diag = 0.0;

    // Retrieve previous value for velocity.
    fe_values[velocity].get_function_values(solution_old, velocity_old_loc);

    for (unsigned int q = 0; q < n_q; ++q)
    {
      // Get the evaluation of the forcing term at the current quadrature point.
      Vector<double> forcing_term_loc(dim);
      forcing_term.vector_value(fe_values.quadrature_point(q),
                                forcing_term_loc);
      Tensor<1, dim> forcing_term_tensor;
      for (unsigned int d = 0; d < dim; ++d)
        forcing_term_tensor[d] = forcing_term_loc[d];

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {
          // Velocity Mass contribution term in the momentum equation.
          cell_matrix(i, j) +=
              fe_values[velocity].value(i, q) *
              fe_values[velocity].value(j, q) /
              (1.0 / (double)settings.get_time_steps_per_second()) * // = deltaT
              fe_values.JxW(q);

          // Stiffness contribution term in the momentum equation.
          cell_matrix(i, j) +=
              get_current_coeff_nu() *
              scalar_product(fe_values[velocity].gradient(i, q),
                             fe_values[velocity].gradient(j, q)) *
              fe_values.JxW(q);

          // Linearized Convection term contribution in the momentum equation.
          cell_matrix(i, j) +=
              fe_values[velocity].value(i, q) *
              (fe_values[velocity].gradient(j, q) *
               velocity_old_loc[q]) *
              fe_values.JxW(q);

          // Pressure term in the momentum equation.
          cell_matrix(i, j) -= fe_values[velocity].divergence(i, q) *
                               fe_values[pressure].value(j, q) *
                               fe_values.JxW(q);

          // Pressure term in the continuity equation.
          cell_matrix(i, j) -= fe_values[velocity].divergence(j, q) *
                               fe_values[pressure].value(i, q) *
                               fe_values.JxW(q);
        }

        // Local contributions to the diagonal of the F matrix:

        // Velocity Mass contribution term in the momentum equation.
        cell_diag(i) +=
            fe_values[velocity].value(i, q) *
            fe_values[velocity].value(i, q) /
            (1.0 / (double)settings.get_time_steps_per_second()) * // = deltaT
            fe_values.JxW(q);

        // Stiffness contribution term in the momentum equation.
        cell_diag(i) +=
            get_current_coeff_nu() *
            scalar_product(fe_values[velocity].gradient(i, q),
                           fe_values[velocity].gradient(i, q)) *
            fe_values.JxW(q);

        // Linearized Convection term contribution in the momentum equation.
        cell_diag(i) +=
            fe_values[velocity].value(i, q) *
            (fe_values[velocity].gradient(i, q) *
             velocity_old_loc[q]) *
            fe_values.JxW(q);

        // Local contributions to the residual vector:

        // Forcing term contribution.
        cell_rhs(i) +=
            scalar_product(forcing_term_tensor,
                           fe_values[velocity].value(i, q)) *
            fe_values.JxW(q);

        // Old-solution term contribution.
        cell_rhs(i) +=
            fe_values[velocity].value(i, q) *
            velocity_old_loc[q] /
            (1.0 / (double)settings.get_time_steps_per_second()) * // deltaT
            fe_values.JxW(q);
      }
    }

    // Boundary integral for Neumann BCs (outlet)
    if (cell->at_boundary())
    {
      for (unsigned int f = 0; f < cell->n_faces(); ++f)
      {
        if (cell->face(f)->at_boundary() &&
            std::find(surfaces_outlets.begin(), surfaces_outlets.end(), cell->face(f)->boundary_id()) != surfaces_outlets.end())
        {
          fe_face_values.reinit(cell, f);

          for (unsigned int q = 0; q < n_q_face; ++q)
          {
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              cell_rhs(i) -=
                  settings.get_outlet_pressure() *
                  scalar_product(fe_face_values.normal_vector(q),
                                 fe_face_values[velocity].value(i,
                                                                q)) *
                  fe_face_values.JxW(q);
            }
          }
        }
      }
    }

    cell->get_dof_indices(dof_indices);

    jacobian_matrix.add(dof_indices, cell_matrix);
    residual_vector.add(dof_indices, cell_rhs);
    D_inv.add(dof_indices, cell_diag);
  }

  jacobian_matrix.compress(VectorOperation::add);
  residual_vector.compress(VectorOperation::add);
  D_inv.compress(VectorOperation::add);

  // Invert D_inv.
  for (unsigned int i = 0; i < D_inv.size(); ++i)
  {
    if (D_inv.in_local_range(i) && D_inv(i) != 0.0)
      D_inv(i) = 1.0 / D_inv(i);
  }

  // Dirichlet boundary conditions.
  {
    std::map<types::global_dof_index, double> boundary_values;
    std::map<types::boundary_id, const Function<dim> *> boundary_functions;

    Functions::ZeroFunction<dim> zero_function(dim + 1);

    // We interpolate first the inlet velocity condition alone, then the wall
    // condition alone, so that the latter "win" over the former where the two
    // boundaries touch.
    // That is, only if the iteration of the Newton method is the first one.
    // Otherwise the Dirichet BCs for the inlet surface in the other steps are all
    // zero (no increase in the inlet velocity).
    for (const unsigned int &surf : surfaces_inlets)
    {
      boundary_functions[surf] = &inlet_velocity;
    }

    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values,
                                             ComponentMask(
                                                 {true, true, true, false}));
    boundary_functions.clear();

    for (const unsigned int &surf : surfaces_walls)
    {
      boundary_functions[surf] = &zero_function;
    }
    // FIXME REMOVEME
    for (const unsigned int &surf : surfaces_free_slip)
    {
      boundary_functions[surf] = &zero_function;
    }

    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values,
                                             ComponentMask(
                                                 {true, true, true, false}));

    MatrixTools::apply_boundary_values(
        boundary_values, jacobian_matrix, solution_owned, residual_vector, false);
  }
}

void NavierStokes::solve_time_step()
{
  pcout << "-----------------------------------------------" << std::endl;
  pcout << "Solving the linear system" << std::endl;

  SolverControl solver_control(settings.get_max_solver_iteration_amt(), settings.get_desired_solver_precision() * residual_vector.l2_norm());

  SolverGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);

  // Custom Defined Preconditioner
  PreconditionSIMPLE preconditioner;
  preconditioner.initialize(settings.get_preconditioner_coeff_alpha(), jacobian_matrix, D_inv.block(0));

  solver.solve(jacobian_matrix, solution_owned, residual_vector, preconditioner);
  pcout << "  " << solver_control.last_step() << " GMRES iterations"
        << std::endl;

  solution = solution_owned;
}

void NavierStokes::solve()
{
  pcout << "===============================================" << std::endl;

  inlet_velocity.set_time_step(current_time_step);
  forcing_term.set_time(0.0);
  initial_solution.set_time(0.0);

  double cur_reynolds_number = estimate_reynolds_number();
  pcout << "Reynold's number estimation at the start: Re = " << cur_reynolds_number << " (" << (cur_reynolds_number <= max_reynolds_number_before_continuation ? "<=" : ">")
        << " " << max_reynolds_number_before_continuation << ")" << std::endl;

  if (cur_reynolds_number > max_reynolds_number_before_continuation)
  {
    // Do continuation method
    pcout << "Continuation method not implemented yet..." << std::endl;
  }
  else
  {
    pcout << "No need for the continuation method." << std::endl;
    // Apply the initial condition.

    pcout << "Applying the initial condition" << std::endl;

    // Here for the initial solution start with the velocity equal to zero.
    VectorTools::interpolate(dof_handler, initial_solution, solution_owned);
    solution = solution_owned;

    // Output the initial solution.
    // output(0, 0.0);
  }

  pcout << "===============================================" << std::endl;

  for (unsigned int t = 0; t < settings.get_time_steps_pre_ramp() + settings.get_time_steps_ramp() + settings.get_time_steps_post_ramp(); ++t)
  {
    pcout << "Current time: " << (double)current_time_step / (double)settings.get_time_steps_per_second()
          << " (time step: " << current_time_step << "/"
          << settings.get_time_steps_pre_ramp() + settings.get_time_steps_ramp() + settings.get_time_steps_post_ramp() << ")" << std::endl;
    pcout << "\tCurrently ";
    if (current_time_step < settings.get_time_steps_pre_ramp())
      pcout << "holding inlet velocity to minimum before ramping up.";
    else if (current_time_step >= settings.get_time_steps_pre_ramp() + settings.get_time_steps_ramp())
      pcout << "holding inlet velocity to maximum after ramping down.";
    else
      pcout << "ramping up inlet velocity.";

    pcout << std::endl;

    Point<dim> temp;
    Vector<double> cur_velocity(dim + 1);
    inlet_velocity.vector_value(temp, cur_velocity);

    pcout << "\tThis time step's velocity: " << cur_velocity[0] << ", "
          << cur_velocity[1] << ", " << cur_velocity[2] << std::endl;
    pcout << "\tThis time step's Reynolds number: " << estimate_reynolds_number() << std::endl;

    solution_old = solution;

    assemble_system();
    solve_time_step();
    output();
    increment_time_step();
  }
}

/*void NavierStokes::solve_newton()
{
  pcout << "===============================================" << std::endl;

  const unsigned int n_max_iters = settings.get_max_newton_iteration_amt();
  const double residual_tolerance = settings.get_desired_newton_precision();

  unsigned int n_iter = 0;
  double residual_norm = residual_tolerance + 1;
  bool initial_step = true;

  while (n_iter < n_max_iters && residual_norm > residual_tolerance)
  {
    assemble_system(initial_step);
    initial_step = false;
    residual_norm = residual_vector.l2_norm();

    // pcout << "Current velocity: " << //get_current_inlet_velocity(0) << ", "
    //       << get_current_inlet_velocity(1) << ", " << get_current_inlet_velocity(2) << std::endl;

    pcout << "Newton iteration " << n_iter << "/" << n_max_iters
          << " - ||r|| = " << std::scientific << std::setprecision(6)
          << residual_norm << std::flush;

    // We actually solve the system only if the residual is larger than the
    // tolerance.
    if (residual_norm > residual_tolerance)
    {
      solve_newton_step();

      solution_owned += delta_owned;
      solution = solution_owned;
    }
    else
    {
      pcout << " < tolerance" << std::endl;
    }

    ++n_iter;
  }

  pcout << "===============================================" << std::endl;
}*/

void NavierStokes::output() const
{
  pcout << "===============================================" << std::endl;

  dealii::Timer output_timer(MPI_COMM_WORLD, true);

  DataOut<dim> data_out;

  // data_out.attach_dof_handler(dof_handler);

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
          dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);
  std::vector<std::string> names = {"velocity",
                                    "velocity",
                                    "velocity",
                                    "pressure"};

  data_out.add_data_vector(dof_handler,
                           solution,
                           names,
                           data_component_interpretation);

  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  std::string temp = std::to_string(settings.get_time_steps_pre_ramp() + settings.get_time_steps_ramp() + settings.get_time_steps_post_ramp());

  const std::string output_file_name = settings.get_out_file_name() + "_" +
                                       std::string(1 + temp.size() - std::to_string(current_time_step).size(), '0') +
                                       std::to_string(current_time_step);

  DataOutBase::DataOutFilter data_filter(
      DataOutBase::DataOutFilterFlags(/*filter_duplicate_vertices = */ false,
                                      /*xdmf_hdf5_output = */ true));
  data_out.write_filtered_data(data_filter);
  data_out.write_hdf5_parallel(data_filter,
                               output_file_name + ".h5",
                               MPI_COMM_WORLD);

  std::vector<XDMFEntry> xdmf_entries({data_out.create_xdmf_entry(
      data_filter, output_file_name + ".h5", (double)current_time_step / (double)settings.get_time_steps_per_second(), MPI_COMM_WORLD)});
  data_out.write_xdmf_file(xdmf_entries,
                           output_file_name + ".xdmf",
                           MPI_COMM_WORLD);

  output_timer.stop();

  pcout << "Output written to " << output_file_name << std::endl;
  pcout << "===============================================" << std::endl;
  if (mpi_rank == 0)
  {
    std::cout.precision(5);
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    pcout << "Output time: " << output_timer.wall_time() << "s" << std::endl;
    pcout << "===============================================" << std::endl;
  }
}