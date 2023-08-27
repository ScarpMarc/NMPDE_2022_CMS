#include "Problem.hpp"
#include "SimulationSettings.hpp"

// Main function.
int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  ns_sim_settings::SimulationSettings sim_settings;
  std::string settings_file_name;

  int mpi_ID;
  int mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_ID);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  if (mpi_ID == 0)
  {
    std::cout << "Running on " << mpi_size << " processes." << std::endl;
  }

  if (argc > 1)
  {
    settings_file_name = argv[1];
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

  //dealii::Timer timer(MPI_COMM_WORLD, true);

  NavierStokes problem(sim_settings);

  try
  {
    problem.setup();
    problem.solve();
    // Output is handled at the end of each time step.
  }
  catch (const std::exception &e)
  {
    if (mpi_ID == 0)
    {
      std::cerr << e.what() << '\n';
    }
  }

  // Prints timer info
  problem.finalize();

  /*timer.stop();
  if (mpi_ID == 0)
  {
    std::cout.precision(5);
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout << "Total runtime: " << timer.wall_time() << "s" << std::endl;
  }*/

  return 0;
}