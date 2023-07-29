#ifndef SIM_SETTINGS_HPP
#define SIM_SETTINGS_HPP

#include <string>
#include <array>
#include <deal.II/base/parameter_handler.h>

// This undefined text serves as a flag in the data structure to indicate that
// the variable should not be loaded from file.
#define INTERNAL_ONLY

namespace ns_sim_settings
{
    class SimulationSettings
    {
    public:
        SimulationSettings(); 

        SimulationSettings(const SimulationSettings &other) : SimulationSettings(other.file_name,
                                                                                 other.degree_velocity,
                                                                                 other.degree_pressure,
                                                                                 other.coeff_nu,
                                                                                 other.coeff_nu_start,
                                                                                 other.coeff_nu_ramp_down_times,
                                                                                 other.coeff_rho,
                                                                                 other.inlet_velocity,
                                                                                 other.inlet_velocity_increment_per_step,
                                                                                 other.outlet_pressure,
                                                                                 other.max_solver_iteration_amt,
                                                                                 other.desired_solver_precision,
                                                                                 other.max_newton_iteration_amt,
                                                                                 other.desired_newton_precision,
                                                                                 other.theta,
                                                                                 other.total_time_steps,
                                                                                 other.time_steps_per_second) {}

        SimulationSettings(const std::string file_name,
                           const unsigned int degree_velocity,
                           const unsigned int degree_pressure,
                           const double coeff_nu,
                           const double coeff_nu_start,
                           const unsigned int coeff_nu_ramp_down_times,
                           const double coeff_rho,
                           const std::array<double, 3> inlet_velocity,
                           const std::array<double, 3> inlet_velocity_increment_per_step,
                           const double outlet_pressure,
                           const unsigned int max_solver_iteration_amt,
                           const double desired_solver_precision,
                           const unsigned int max_newton_iteration_amt,
                           const double desired_newton_precision,
                           const double theta,
                           const unsigned long total_time_steps,
                           const unsigned long time_steps_per_second);

        void read_data(const std::string &filename);
        void print() const;

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
        double advised_ramp_down_step_size;

        dealii::ParameterHandler prm;
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