#include "Problem.hpp"
#include "SimulationSettings.hpp"

// Main function.
int main(int argc, char *argv[])
{
  ns_sym_settings::SimulationSettings sim_settings = ns_sym_settings::base_sim_settings;

  if (argc > 1)
  {
    // Load from file in the future. For now throw error
    // sym_settings = load_settings(argv[1]);
    std::cout << "Error! File loading will be supported in the future" << std::endl;
    return -1;
  }

  ns_sym_settings::print_simulation_settings(sim_settings);

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  Problem problem(sym_settings);

  problem.setup();
  problem.assemble();
  problem.solve();
  problem.output();

  return 0;
}