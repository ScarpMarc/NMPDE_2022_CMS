#include "Problem.hpp"
#include "SimulationSettings.hpp"

// Main function.
int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  ns_sim_settings::SimulationSettings sim_settings;
  std::string settings_file_name;

  settings_file_name = argv[1];

  int mpi_ID;
  int mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_ID);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  if (argc > 1)
  {
    if (mpi_ID == 0)
    {
      std::cout << "Loading settings from file: " << settings_file_name << std::endl;
    }
    sim_settings.read_data(settings_file_name);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (mpi_ID == 0)
  {
    sim_settings.print();
  }

  NavierStokes problem(sim_settings);

  problem.setup();
  problem.solve_newton();
  problem.output();

  return 0;
}