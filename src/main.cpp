#include "Problem.hpp"
#include "SimulationSettings.hpp"

SimulationSettings base_sym_settings = {
    "Ptero_reducedx2.msh", // std::string file_name;
    2,                          // unsigned int degree_velocity;
    1,                          // unsigned int degree_pressure;
    1.0,                        // double coeff_nu;
    {0.0, 10.0, 0.0},           // std::array<double, 3> inlet_velocity;
    10.0,                       // double outlet_pressure;
    5000,                       // unsigned int max_solver_iteration_amt;
    1e-6,                       // double desired_solver_precision;
    0.5                         // double theta;
};

// Main function.
int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  SimulationSettings sym_settings = base_sym_settings;

  if (argc > 1)
  {
    // Load from file in the future. For now throw error
    // sym_settings = load_settings(argv[1]);
    std::cout << "Error! File loading will be supported in the future" << std::endl;
    return -1;
  }

  // Output settings
  std::cout << "File name: " << sym_settings.file_name << std::endl; //"half_ptero_reducedx2.msh",//std::string file_name;
  std::cout << std::endl;
  
  std::cout << "nu coefficient: " << sym_settings.coeff_nu << std::endl;              // 1.0,//double coeff_nu;
  std::cout << "Outlet pressure (Pa): " << sym_settings.outlet_pressure << std::endl; // 10.0,//double outlet_pressure;
  std::cout << "Inlet velocity components (m/s): ";
  for (unsigned int i = 0; i < 3; ++i)
    std::cout << sym_settings.inlet_velocity[i] << " "; //{0.0, 10.0, 0.0},//std::array<double, 3> inlet_velocity;
  std::cout << std::endl;

  std::cout << std::endl;

  std::cout << "Velocity degree: " << sym_settings.degree_velocity << std::endl; // 2,//unsigned int degree_velocity;
  std::cout << "Pressure degree: " << sym_settings.degree_pressure << std::endl; // 1,//unsigned int degree_pressure;
  std::cout << std::endl;

  std::cout << "Max. solver iterations: " << sym_settings.max_solver_iteration_amt << std::endl;   // 5000,//unsigned int max_solver_iteration_amt;
  std::cout << "Desired solver precision: " << sym_settings.desired_solver_precision << std::endl; // 1e-6,//double desired_solver_precision;
  std::cout << "Theta for Theta Method: " << sym_settings.theta << std::endl;                      // 0.5//double theta;

  Problem problem(sym_settings);

  problem.setup();
  problem.assemble();
  problem.solve();
  problem.output();

  return 0;
}