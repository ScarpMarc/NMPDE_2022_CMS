#include "SimulationSettings.hpp"
#include <iostream>

namespace ns_sym_settings
{
    void print_simulation_settings(const SimulationSettings &sim_settings)
    {
        // Output settings
        std::cout << "File name: " << sim_settings.file_name << std::endl; //"half_ptero_reducedx2.msh",//std::string file_name;
        std::cout << std::endl;

        std::cout << "nu coefficient: " << sim_settings.coeff_nu << std::endl;              // 1.0,//double coeff_nu;
        std::cout << "Outlet pressure (Pa): " << sim_settings.outlet_pressure << std::endl; // 10.0,//double outlet_pressure;
        std::cout << "Inlet velocity components (m/s): ";
        for (unsigned int i = 0; i < 3; ++i)
            std::cout << sim_settings.inlet_velocity[i] << " "; //{0.0, 10.0, 0.0},//std::array<double, 3> inlet_velocity;
        std::cout << std::endl;

        std::cout << std::endl;

        std::cout << "Velocity degree: " << sim_settings.degree_velocity << std::endl; // 2,//unsigned int degree_velocity;
        std::cout << "Pressure degree: " << sim_settings.degree_pressure << std::endl; // 1,//unsigned int degree_pressure;
        std::cout << std::endl;

        std::cout << "Max. solver iterations: " << sim_settings.max_solver_iteration_amt << std::endl;   // 5000,//unsigned int max_solver_iteration_amt;
        std::cout << "Desired solver precision: " << sim_settings.desired_solver_precision << std::endl; // 1e-6,//double desired_solver_precision;
        std::cout << "Theta for Theta Method: " << sim_settings.theta << std::endl;                      // 0.5//double theta;
    }
}