#ifndef SIM_SETTINGS_HPP
#define SIM_SETTINGS_HPP

#include <string>
#include <array>

// This undefined text serves as a flag in the data structure to indicate that
// the variable should not be loaded from file.
#define INTERNAL_ONLY

namespace ns_sim_settings
{
    class SimulationSettings
    {
    public:
        void read_data(const std::string &filename);
        void print();

        std::string file_name;
        unsigned int degree_velocity;
        unsigned int degree_pressure;

        double coeff_nu;                       // Dynamic viscosity [Pa*s]
        double coeff_nu_start;                 // Start of dynamic viscosity for iterative "ramp-down" process [Pa*s]
        unsigned int coeff_nu_ramp_down_times; // Amount of times to ramp down viscosity

        double coeff_rho; // Density [kg/m^3]

        std::array<double, 3> inlet_velocity;
        std::array<double, 3> inlet_velocity_increment_per_step;
        double outlet_pressure;

        unsigned int max_solver_iteration_amt;
        double desired_solver_precision;

        unsigned int max_newton_iteration_amt;
        double desired_newton_precision;

        double theta;

        unsigned long total_time_steps;
        unsigned long time_steps_per_second;

    protected:
    // Variables here are determined automatically and/or serve just as info for the user.
        double ramp_down_step_size;
        const double advised_ramp_down_step_size;

        ParameterHandler prm;
    };

    /*const SimulationSettings base_sim_settings = {
        "Ptero_full.msh", // std::string file_name;
        2,                // unsigned int degree_velocity;
        1,                // unsigned int degree_pressure;
        1.825e-5,         // double coeff_nu, Dynamic viscosity
        2.5e-3,           // double coeff_nu_start; // Start of dynamic viscosity for iterative "ramp-down" process
        10,               // unsigned int ramp_down_times; // Number of times to ramp down viscosity;
        {0.0, 10.0, 0.0}, // std::array<double, 3> inlet_velocity;
        10.0,             // double outlet_pressure;
        5000,             // unsigned int max_solver_iteration_amt;
        1e-6,             // double desired_solver_precision;
        1000,             // unsigned int max_newton_iteration_amt;
        1e-6,             // double desired_newton_precision;
        0.5               // double theta;
    };*/
}

#endif