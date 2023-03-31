#ifndef SIM_SETTINGS_HPP
#define SIM_SETTINGS_HPP

#include <string>
#include <array>

namespace ns_sim_settings
{
    struct SimulationSettings
    {
        std::string file_name;
        unsigned int degree_velocity;
        unsigned int degree_pressure;

        double coeff_nu;
        double gamma;
        std::array<double, 3> inlet_velocity;
        double outlet_pressure;

        unsigned int max_solver_iteration_amt;
        double desired_solver_precision;

        unsigned int max_newton_iteration_amt;
        double desired_newton_precision;

        double theta;

        unsigned long total_time_steps;
        unsigned long time_steps_per_second;
    };

    const SimulationSettings base_sim_settings = {
        "Ptero_full.msh", // std::string file_name;
        2,                      // unsigned int degree_velocity;
        1,                      // unsigned int degree_pressure;
        1.0,                    // double coeff_nu;
        1.0,                    // double gamma;
        {0.0, 10.0, 0.0},       // std::array<double, 3> inlet_velocity;
        10.0,                   // double outlet_pressure;
        5000,                   // unsigned int max_solver_iteration_amt;
        1e-6,                   // double desired_solver_precision;
        1000,                   //unsigned int max_newton_iteration_amt;
        1e-6,                   //double desired_newton_precision;
        0.5,                     // double theta;
        10,                     //unsigned long total_time_steps;
        10                      //unsigned long time_steps_per_second;
    };

    void print_simulation_settings(const SimulationSettings &sim_settings);
}

#endif