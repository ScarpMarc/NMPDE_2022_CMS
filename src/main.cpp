#include "Problem.hpp"
#include "SimulationSettings.hpp"

// Main function.
int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  ns_sim_settings::SimulationSettings sim_settings;

  int mpi_ID;

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_ID);
  if (mpi_ID == 0)
  {
    if (argc > 1)
    {
      // Load from file in the future. For now throw an error
      // sym_settings = load_settings(argv[1]);
      std::cout << "Error! File loading will be supported in the future" << std::endl;
      return -1;
    }

    sim_settings.print();
  }

  Problem problem(sim_settings);

  problem.setup();
  problem.assemble();
  problem.solveNewtonMethod();
  //problem.solve();
  problem.output();

  return 0;
}