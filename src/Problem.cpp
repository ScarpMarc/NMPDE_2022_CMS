#include "Problem.hpp"
#include "PreconditionBlockTriangular.hpp"

void Problem::setup()
{
  // Create the mesh.
  {
    pcout << "Initializing the mesh" << std::endl;

    Triangulation<dim> mesh_serial;

    GridIn<dim> grid_in;
    grid_in.attach_triangulation(mesh_serial);

    /*const std::string mesh_file_name =
      "../mesh/mesh-step-" + std::to_string(N + 1) + ".msh";*/

    std::ifstream grid_in_file(simulation_settings.file_name);
    grid_in.read_msh(grid_in_file);

    GridTools::partition_triangulation(mpi_size, mesh_serial);
    const auto construction_data = TriangulationDescription::Utilities::
        create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
    mesh.create_triangulation(construction_data);

    pcout << "  Number of elements = " << mesh.n_global_active_cells()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space.
  {
    pcout << "Initializing the finite element space" << std::endl;

    const FE_SimplexP<dim> fe_scalar_velocity(simulation_settings.degree_velocity);
    const FE_SimplexP<dim> fe_scalar_pressure(simulation_settings.degree_pressure);
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
    for (unsigned int c = 0; c < dim + 1; ++c)
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
    sparsity_pressure_mass.compress();

    pcout << "  Initializing the matrices" << std::endl;
    jacobian_matrix.reinit(sparsity);
    pressure_mass.reinit(sparsity_pressure_mass);

    pcout << "  Initializing the system right-hand side" << std::endl;
    residual_vector.reinit(block_owned_dofs, MPI_COMM_WORLD);
    pcout << "  Initializing the solution vector" << std::endl;
    solution_owned.reinit(block_owned_dofs, MPI_COMM_WORLD);
    solution.reinit(block_owned_dofs, block_relevant_dofs, MPI_COMM_WORLD);
    delta_owned.reinit(block_owned_dofs, MPI_COMM_WORLD);
  }
}

/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/

void Problem::assemble()
{
  pcout << "===============================================" << std::endl;
  pcout << "Assembling the system" << std::endl;

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
  FullMatrix<double> cell_pressure_mass_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  jacobian_matrix = 0.0;
  residual_vector = 0.0;
  // This one just for the construction of the preconditioner
  // (Should be built at every step of the Newton method)
  pressure_mass = 0.0;

  // We use these vectors to store the old solution (i.e. at previous Newton
  // iteration) and its gradient on quadrature nodes of the current cell:

  // Tensor of rank 1 to store the local velocity solution at previous Newton step
  //(vectorial field)
  std::vector<Tensor<1, dim>> velocity_loc(n_q);
  // Tensor of rank 2 to store the local velocity solution at previous Newton step
  //(the grad of a vectorial field -> jacobian)
  std::vector<Tensor<2, dim>> velocity_gradient_loc(n_q);

  // Double to store the local pressure solution at previous Newton step
  std::vector<double> pressure_loc(n_q);

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (!cell->is_locally_owned())
      continue;

    fe_values.reinit(cell);

    cell_matrix = 0.0;
    cell_rhs = 0.0;
    cell_pressure_mass_matrix = 0.0;

    // We need to compute the Jacobian matrix and the residual for current
    // cell. This requires knowing the value and the gradient of u^{(k)}
    // (stored inside solution) on the quadrature nodes of the current
    // cell. This can be accomplished through
    // FEValues::get_function_values and FEValues::get_function_gradients.
    fe_values[velocity].get_function_values(solution, velocity_loc);
    fe_values[velocity].get_function_gradients(solution, velocity_gradient_loc);
    fe_values[pressure].get_function_values(solution, pressure_loc);

    for (unsigned int q = 0; q < n_q; ++q)
    {
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
          // Viscosity term.
          cell_matrix(i, j) +=
              simulation_settings.coeff_nu *
              scalar_product(fe_values[velocity].gradient(j, q),
                             fe_values[velocity].gradient(i, q)) *
              fe_values.JxW(q);

          // Convective term
          cell_matrix(i, j) +=
              (velocity_gradient_loc[q] *
                   fe_values[velocity].value(j, q) *
                   fe_values[velocity].value(i, q) +
               fe_values[velocity].gradient(j, q) *
                   velocity_loc[q] *
                   fe_values[velocity].value(i, q)) *
              fe_values.JxW(q);

          // Pressure term in the momentum equation.
          cell_matrix(i, j) -= fe_values[velocity].divergence(j, q) *
                               fe_values[pressure].value(i, q) *
                               fe_values.JxW(q);

          // Pressure term in the continuity equation.
          cell_matrix(i, j) -= fe_values[pressure].value(j, q) *
                               fe_values[velocity].divergence(i, q) *
                               fe_values.JxW(q);

          // Term given by the special-case preconditioer
          cell_matrix(i, j) -= gamma *
                               fe_values[velocity].divergence(j, q) *
                               fe_values[velocity].divergence(i, q) *
                               fe_values.JxW(q);

          // Pressure mass matrix.
          cell_pressure_mass_matrix(i, j) +=
              fe_values[pressure].value(i, q) *
              fe_values[pressure].value(j, q) / simulation_settings.coeff_nu * fe_values.JxW(q);
        }

        // Get the divergence of the previous

        double velocity_divergence_loc =
            trace(velocity_gradient_loc[q]);

        // Forcing term.

        cell_rhs(i) +=
            ( // Viscosity term
                simulation_settings.coeff_nu *
                    scalar_product(velocity_gradient_loc[q],
                                   fe_values[velocity].gradient(i, q)) *
                    // Convolution term
                    -velocity_gradient_loc[q] *
                    velocity_loc[q] *
                    fe_values[velocity].value(i, q)
                // Pressure term
                +
                pressure_loc[q] *
                    fe_values[velocity].divergence(i, q) +
                velocity_divergence_loc *
                    fe_values[pressure].value(i, q) -
                gamma *
                    velocity_divergence_loc *
                    fe_values[velocity].divergence(i, q) -
                scalar_product(forcing_term_tensor,
                               fe_values[velocity].value(i, q))) *
            fe_values.JxW(q);
      }
    }

    // Boundary integral for Neumann BCs.
    if (cell->at_boundary())
    {
      for (unsigned int f = 0; f < cell->n_faces(); ++f)
      {
        if (cell->face(f)->at_boundary() &&
            cell->face(f)->boundary_id() == 2)
        {
          fe_face_values.reinit(cell, f);

          for (unsigned int q = 0; q < n_q_face; ++q)
          {
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              cell_rhs(i) +=
                  -simulation_settings.outlet_pressure *
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
    pressure_mass.add(dof_indices, cell_pressure_mass_matrix);
  }

  jacobian_matrix.compress(VectorOperation::add);
  residual_vector.compress(VectorOperation::add);
  pressure_mass.compress(VectorOperation::add);

  // Dirichlet boundary conditions.
  {
    std::map<types::global_dof_index, double> boundary_values;
    std::map<types::boundary_id, const Function<dim> *> boundary_functions;

    // We interpolate first the inlet velocity condition alone, then the wall
    // condition alone, so that the latter "win" over the former where the two
    // boundaries touch.
    boundary_functions[1] = &inlet_velocity;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values,
                                             ComponentMask(
                                                 {true, true, true, false}));

    boundary_functions.clear();
    Functions::ZeroFunction<dim> zero_function(dim + 1);
    boundary_functions[0] = &zero_function; // Body
    boundary_functions[3] = &zero_function; // Slip
    boundary_functions[4] = &zero_function; // Wall

    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values,
                                             ComponentMask(
                                                 {true, true, true, false}));

    MatrixTools::apply_boundary_values(
        boundary_values, jacobian_matrix, solution, residual_vector, false);
  }
}

/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/

void Problem::solveLinearSystem()
{
  /*
  // IN : jacobian_matrix, residual_vector, preconditionr
  // OUT: delta_owned

  // inizializzo il controllo del solver
  SolverControl solver_control(1000, 1e-6 * residual_vector.l2_norm());

  // inizializzo il solver
  SolverCG<TrilinosWrappers::MPI::BlockVector> solver(solver_control);

  // definisco un precondizionatore
  TrilinosWrappers::PreconditionSSOR      preconditioner;
  preconditioner.initialize(
    jacobian_matrix, TrilinosWrappers::PreconditionSSOR::AdditionalData(1.0));

  solver.solve(jacobian_matrix,
               delta_owned,
               residual_vector,
               preconditioner);
  pcout << "  " << solver_control.last_step() << " CG iterations" << std::endl;
  */

  pcout << "===============================================" << std::endl;

  SolverControl solver_control(simulation_settings.max_solver_iteration_amt, simulation_settings.desired_solver_precision * residual_vector.l2_norm());

  SolverGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);

  // PreconditionBlockDiagonal preconditioner;
  // preconditioner.initialize(system_matrix.block(0, 0),
  //                           pressure_mass.block(1, 1));

  PreconditionBlockTriangular preconditioner;
  preconditioner.initialize(jacobian_matrix.block(0, 0),
                            pressure_mass.block(1, 1),
                            jacobian_matrix.block(0, 1),
                            simulation_settings.coeff_nu,
                            gamma);

  pcout << "Solving the linear system" << std::endl;
  solver.solve(jacobian_matrix,
               delta_owned,
               residual_vector,
               preconditioner);
  pcout << "  " << solver_control.last_step() << " GMRES iterations"
        << std::endl;
}

void Problem::solveNewtonMethod()
{
  pcout << "===============================================" << std::endl;
  // parametri
  const unsigned int n_max_iters = simulation_settings.max_newton_iteration_amt;
  const double residual_tolerance = simulation_settings.desired_newton_precision;

  unsigned int n_iter = 0;
  double residual_norm = residual_tolerance + 1;

  // We apply the boundary conditions to the initial guess (which is stored in
  // solution_owned and solution).
  {
    IndexSet dirichlet_dofs = DoFTools::extract_boundary_dofs(dof_handler);

    // function_g.set_time(time);

    TrilinosWrappers::MPI::BlockVector vector_dirichlet(solution_owned);
    VectorTools::interpolate(dof_handler,
                             inlet_velocity,
                             vector_dirichlet);

    for (const auto &idx : dirichlet_dofs)
      solution_owned[idx] = vector_dirichlet[idx];

    solution_owned.compress(VectorOperation::insert);
    solution = solution_owned;
  }

  while (n_iter < n_max_iters && residual_norm > residual_tolerance)
  {
    pcout << "Newton: back to assembly" << std::endl;
    assemble();                                // riassemblo
    residual_norm = residual_vector.l2_norm(); // calcolo la norma residuo

    pcout << "  Newton iteration " << n_iter << "/" << n_max_iters
          << " - ||r|| = " << std::scientific << std::setprecision(6)
          << residual_norm << std::flush << std::endl;

    // We actually solve the system only if the residual is larger than the
    // tolerance.
    if (residual_norm > residual_tolerance)
    {
      solveLinearSystem();

      solution_owned += delta_owned;
      solution = solution_owned;
    }
    else
    {
      pcout << " < tolerance" << std::endl;
    }

    ++n_iter;
  }
}

void Problem::output()
{
  pcout << "===============================================" << std::endl;

  DataOut<dim> data_out;

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

  const std::string output_file_name = "output-Ptero";

  DataOutBase::DataOutFilter data_filter(
      DataOutBase::DataOutFilterFlags(/*filter_duplicate_vertices = */ false,
                                      /*xdmf_hdf5_output = */ true));
  data_out.write_filtered_data(data_filter);
  data_out.write_hdf5_parallel(data_filter,
                               output_file_name + ".h5",
                               MPI_COMM_WORLD);

  std::vector<XDMFEntry> xdmf_entries({data_out.create_xdmf_entry(
      data_filter, output_file_name + ".h5", 0, MPI_COMM_WORLD)});
  data_out.write_xdmf_file(xdmf_entries,
                           output_file_name + ".xdmf",
                           MPI_COMM_WORLD);

  pcout << "Output written to " << output_file_name << std::endl;
  pcout << "===============================================" << std::endl;
}