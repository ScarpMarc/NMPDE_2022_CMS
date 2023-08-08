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

        const std::string file_name() const { return prm.get("file_name"); }
        const std::string out_file_name() const { return prm.get("out_file_name"); };
        unsigned int degree_velocity()
        {
            unsigned int output;
            prm.enter_subsection("Finite-element settings");
            {
                output = prm.get_integer("degree_velocity");
            }
            prm.leave_subsection();
            return output;
        }

        unsigned int degree_pressure()
        {
            unsigned int output;
            prm.enter_subsection("Finite-element settings");
            {
                output = prm.get_integer("degree_pressure");
            }
            prm.leave_subsection();
            return output;
        }

        double coeff_nu() // Dynamic viscosity [Pa*s]
        {
            double output;
            prm.enter_subsection("Physical coefficients");
            {
                prm.enter_subsection("Viscosity");
                {
                    output = prm.get_double("coeff_nu");
                }
                prm.leave_subsection();
            }
            prm.leave_subsection();
            return output;
        }
        double coeff_nu_start() // Start of dynamic viscosity for iterative "ramp-down" process [Pa*s]
        {
            double output;
            prm.enter_subsection("Physical coefficients");
            {
                prm.enter_subsection("Viscosity");
                {
                    output = prm.get_double("coeff_nu_start");
                }
                prm.leave_subsection();
            }
            prm.leave_subsection();
            return output;
        }
        unsigned int coeff_nu_ramp_down_times() // Amount of times to ramp down viscosity
        {
            unsigned int output;
            prm.enter_subsection("Physical coefficients");
            {
                prm.enter_subsection("Viscosity");
                {
                    output = prm.get_integer("coeff_nu_ramp_down_times");
                }
                prm.leave_subsection();
            }
            prm.leave_subsection();
            return output;
        }

        double coeff_rho() // Density [kg/m^3]
        {
            double output;
            prm.enter_subsection("Physical coefficients");
            {
                output = prm.get_double("rho");
            }
            prm.leave_subsection();
            return output;
        }

        const std::array<double, 3> inlet_velocity_start()
        {
            std::array<double, 3> output;
            prm.enter_subsection("Inlet velocity");
            {
                prm.enter_subsection("Start");
                {
                    output[0] = prm.get_double("x");
                    output[1] = prm.get_double("y");
                    output[2] = prm.get_double("z");
                }
                prm.leave_subsection();
            }
            prm.leave_subsection();

            return output;
        }

        const std::array<double, 3> inlet_velocity_end()
        {
            std::array<double, 3> output;

            prm.enter_subsection("Inlet velocity");
            {
                prm.enter_subsection("End");
                {
                    output[0] = prm.get_double("x");
                    output[1] = prm.get_double("y");
                    output[2] = prm.get_double("z");
                }
                prm.leave_subsection();
            }
            prm.leave_subsection();

            return output;
        }
        double outlet_pressure() const
        {
            return prm.get_double("outlet_pressure");
        }

        unsigned int max_solver_iteration_amt()
        {
            unsigned int output;
            prm.enter_subsection("Solver settings");
            {
                output = prm.get_integer("max_solver_iteration_amt");
            }
            prm.leave_subsection();
            return output;
        }

        double desired_solver_precision()
        {
            double output;
            prm.enter_subsection("Solver settings");
            {
                output = prm.get_double("desired_solver_precision");
            }
            prm.leave_subsection();
            return output;
        }

        unsigned int max_newton_iteration_amt()
        {
            unsigned int output;
            prm.enter_subsection("Solver settings");
            {
                output = prm.get_integer("max_newton_iteration_amt");
            }
            prm.leave_subsection();
            return output;
        }
        double desired_newton_precision()
        {
            double output;
            prm.enter_subsection("Solver settings");
            {
                output = prm.get_double("desired_newton_precision");
            }
            prm.leave_subsection();
            return output;
        }

        double theta()
        {
            double output;
            prm.enter_subsection("Solver settings");
            {
                output = prm.get_double("theta");
            }
            prm.leave_subsection();
            return output;
        }
        double coeff_relax_gamma()
        {
            double output;
            prm.enter_subsection("Solver settings");
            {
                output = prm.get_double("coeff_relax_gamma");
            }
            prm.leave_subsection();
            return output;
        }

        unsigned long total_time_steps()
        {
            unsigned long output;
            prm.enter_subsection("Time settings");
            {
                output = prm.get_integer("total_time_steps");
            }
            prm.leave_subsection();
            return output;
        }
        unsigned long time_steps_per_second()
        {
            unsigned long output;
            prm.enter_subsection("Time settings");
            {
                output = prm.get_integer("time_steps_per_second");
            }
            prm.leave_subsection();
            return output;
        }

    protected:
        // Variables here are determined automatically and/or serve just as info for the user.
        double ramp_down_step_size;
        double advised_ramp_down_step_size;

        void initialise_prm();

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