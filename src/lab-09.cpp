#include "Stokes.hpp"

std::string mesh_file_name_glob;

// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  mesh_file_name_glob = argv[1];

  const unsigned int N               = 8;
  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;

  Stokes problem(N, degree_velocity, degree_pressure);

  problem.setup();
  problem.assemble();
  problem.solve();
  problem.output();

  return 0;
}