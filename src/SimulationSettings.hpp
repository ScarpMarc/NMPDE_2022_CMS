#ifndef SIM_SETTINGS_HPP
#define SIM_SETTINGS_HPP

#include <string>
#include <array>

struct SimulationSettings
{
    std::string file_name;
    unsigned int degree_velocity;
    unsigned int degree_pressure;

    double coeff_nu;
    std::array<double, 3> inlet_velocity;
    double outlet_pressure;

    unsigned int max_solver_iteration_amt;
    double desired_solver_precision;

    double theta;
};

#endif  