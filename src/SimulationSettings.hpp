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
        SimulationSettings() { initialise_prm(); }

        SimulationSettings(const std::string &file_name)
        {
            initialise_prm();
            read_data(file_name);
        }

        /* SimulationSettings(const SimulationSettings &other) : SimulationSettings(other.file_name,
                                                                                  other.out_file_name,
                                                                                  other.degree_velocity,
                                                                                  other.degree_pressure,
                                                                                  other.coeff_nu,
                                                                                  other.coeff_nu_start,
                                                                                  other.coeff_nu_ramp_down_times,
                                                                                  other.coeff_rho,
                                                                                  other.inlet_velocity_start,
                                                                                  other.inlet_velocity_end,
                                                                                  other.outlet_pressure,
                                                                                  other.max_solver_iteration_amt,
                                                                                  other.desired_solver_precision,
                                                                                  other.max_newton_iteration_amt,
                                                                                  other.desired_newton_precision,
                                                                                  other.theta,
                                                                                  other.coeff_relax_gamma,
                                                                                  other.total_time_steps,
                                                                                  other.time_steps_per_second) {}*/
        // SimulationSettings(const SimulationSettings &other) : prm(other.prm) {}

        /*SimulationSettings(const std::string file_name,
                           const std::string out_file_name,
                           const unsigned int degree_velocity,
                           const unsigned int degree_pressure,
                           const double coeff_nu,
                           const double coeff_nu_start,
                           const unsigned int coeff_nu_ramp_down_times,
                           const double coeff_rho,
                           const std::array<double, 3> inlet_velocity_start,
                           const std::array<double, 3> inlet_velocity_end,
                           // const std::array<double, 3> inlet_velocity_increment_per_step,
                           const double outlet_pressure,
                           const unsigned int max_solver_iteration_amt,
                           const double desired_solver_precision,
                           const unsigned int max_newton_iteration_amt,
                           const double desired_newton_precision,
                           const double theta,
                           const double coeff_relax_gamma,
                           const unsigned long total_time_steps,
                           const unsigned long time_steps_per_second);*/

        void read_data(const std::string &filename);
        void print() const;

        std::string get_file_name() const { return file_name; }

        std::string get_out_file_name() const { return out_file_name; }

        constexpr unsigned int get_degree_velocity() const { return degree_velocity; }

        constexpr unsigned int get_degree_pressure() const { return degree_pressure; }

        constexpr double get_coeff_nu() const { return coeff_nu; }

        constexpr double get_coeff_nu_start() const { return coeff_nu_start; }

        constexpr unsigned int get_coeff_nu_ramp_down_times() const { return coeff_nu_ramp_down_times; }

        constexpr double get_coeff_rho() const { return coeff_rho; }

        constexpr std::array<double, 3> get_inlet_velocity_start() const { return inlet_velocity_start; }

        constexpr std::array<double, 3> get_inlet_velocity_end() const { return inlet_velocity_end; }

        constexpr double get_outlet_pressure() const { return outlet_pressure; }

        constexpr unsigned int get_max_solver_iteration_amt() const { return max_solver_iteration_amt; }

        constexpr double get_desired_solver_precision() const { return desired_solver_precision; }

        constexpr unsigned int get_max_newton_iteration_amt() const { return max_newton_iteration_amt; }

        constexpr double get_desired_newton_precision() const { return desired_newton_precision; }

        //constexpr double get_theta() const { return theta; }

        constexpr double get_preconditioner_coeff_alpha() const { return preconditioner_coeff_alpha; }

        //constexpr double get_coeff_relax_gamma() const { return coeff_relax_gamma; }

        constexpr unsigned long get_time_steps_pre_ramp() const { return time_steps_pre_ramp; }

        constexpr unsigned long get_time_steps_ramp() const { return time_steps_ramp; }

        constexpr unsigned long get_time_steps_post_ramp() const { return time_steps_post_ramp; }

        constexpr unsigned long get_time_steps_per_second() const { return time_steps_per_second; }

        std::vector<unsigned int> get_surfaces_walls() const { return surfaces_walls; }

        std::vector<unsigned int> get_surfaces_inlets() const { return surfaces_inlets; }

        std::vector<unsigned int> get_surfaces_outlets() const { return surfaces_outlets; }

        std::vector<unsigned int> get_surfaces_free_slip() const { return surfaces_free_slip; }

        double get_characteristic_length() const { return characteristic_length; }

    protected:
        // These two are determined automatically and/or serve just as info for the user.
        double ramp_down_step_size;
        double advised_ramp_down_step_size;

        void initialise_prm();
        void update_variables();

        dealii::ParameterHandler prm;

        std::string file_name;
        std::string out_file_name;
        unsigned int degree_velocity;
        unsigned int degree_pressure;
        double coeff_nu;
        double coeff_nu_start;
        unsigned int coeff_nu_ramp_down_times;
        double coeff_rho;
        std::array<double, 3> inlet_velocity_start;
        std::array<double, 3> inlet_velocity_end;
        double outlet_pressure;
        unsigned int max_solver_iteration_amt;
        double desired_solver_precision;
        unsigned int max_newton_iteration_amt;
        double desired_newton_precision;
        //double theta;
        double preconditioner_coeff_alpha;
        //double coeff_relax_gamma;
        unsigned long time_steps_pre_ramp;
        unsigned long time_steps_ramp;
        unsigned long time_steps_post_ramp;
        unsigned long time_steps_per_second;
        std::vector<unsigned int> surfaces_walls;
        std::vector<unsigned int> surfaces_inlets;
        std::vector<unsigned int> surfaces_outlets;
        std::vector<unsigned int> surfaces_free_slip;
        double characteristic_length;
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