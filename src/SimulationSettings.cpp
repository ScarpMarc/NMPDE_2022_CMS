#include "SimulationSettings.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace dealii;

namespace ns_sim_settings
{
    /*SimulationSettings::SimulationSettings()
        : SimulationSettings(
              "Ptero_full.msh", // file_name("Ptero_full.msh"),
              "Ptero_out",      // file_name("Ptero_full.msh"),
              2,                //        degree_velocity(2),
              1,                //        degree_pressure(1),
              1.825e-5,         //        coeff_nu(1.825e-5),
              2.5e-3,           //        coeff_nu_start(2.5e-3),
              10,               //        coeff_nu_ramp_down_times(10),
              1.204,            //        coeff_rho(1.204),
              {10.0, 0.0, 0.0}, //        inlet_velocity_start({0.0, 10.0, 0.0}),
              {10.0, 0.0, 0.0}, //        inlet_velocity_end({0.0, .1, 0.0}),
              10.0,             //        outlet_pressure(10.0),
              5000,             //        max_solver_iteration_amt(5000),
              1e-6,             //        desired_solver_precision(1e-6),
              1000,             //        max_newton_iteration_amt(1000),
              1e-6,             //        desired_newton_precision(1e-6),
              0.5,              //        theta(0.5),
              1.0,              //        coeff_relax_gamma(1.0),
              100,              //        total_time_steps(100),
              10                //        time_steps_per_second(10)
          )
    {
    }*/

    /*SimulationSettings::SimulationSettings(const std::string file_name,
                                           const std::string out_file_name,
                                           const unsigned int degree_velocity,
                                           const unsigned int degree_pressure,
                                           const double coeff_nu,
                                           const double coeff_nu_start,
                                           const unsigned int coeff_nu_ramp_down_times,
                                           const double coeff_rho,
                                           const std::array<double, 3> inlet_velocity_start,
                                           const std::array<double, 3> inlet_velocity_end,
                                           const double outlet_pressure,
                                           const unsigned int max_solver_iteration_amt,
                                           const double desired_solver_precision,
                                           const unsigned int max_newton_iteration_amt,
                                           const double desired_newton_precision,
                                           const double theta,
                                           const double coeff_relax_gamma,
                                           const unsigned long total_time_steps,
                                           const unsigned long time_steps_per_second) : file_name(file_name),
                                                                                        out_file_name(out_file_name),
                                                                                        degree_velocity(degree_velocity),
                                                                                        degree_pressure(degree_pressure),
                                                                                        coeff_nu(coeff_nu),
                                                                                        coeff_nu_start(coeff_nu_start),
                                                                                        coeff_nu_ramp_down_times(coeff_nu_ramp_down_times),
                                                                                        coeff_rho(coeff_rho),
                                                                                        inlet_velocity_start(inlet_velocity_start),
                                                                                        inlet_velocity_end(inlet_velocity_end),
                                                                                        outlet_pressure(outlet_pressure),
                                                                                        max_solver_iteration_amt(max_solver_iteration_amt),
                                                                                        desired_solver_precision(desired_solver_precision),
                                                                                        max_newton_iteration_amt(max_newton_iteration_amt),
                                                                                        desired_newton_precision(desired_newton_precision),
                                                                                        theta(theta),
                                                                                        coeff_relax_gamma(coeff_relax_gamma),
                                                                                        total_time_steps(total_time_steps),
                                                                                        time_steps_per_second(time_steps_per_second),
                                                                                        ramp_down_step_size(0.0),
                                                                                        advised_ramp_down_step_size(0.0)
    {
        initialise_prm();
    }*/

    void SimulationSettings::initialise_prm()
    {
        prm.declare_entry("file_name",
                          "Ptero_full.msh",
                          Patterns::FileName(),
                          " Mesh file name.");

        prm.declare_entry("out_file_name",
                          "Ptero_out",
                          Patterns::FileName(),
                          " Output file name.");

        prm.enter_subsection("Finite-element settings");
        {
            prm.declare_entry("degree_velocity",
                              "2",
                              Patterns::Integer(1, 3),
                              " Degree of the velocity finite element space. ");
            prm.declare_entry("degree_pressure",
                              "1",
                              Patterns::Integer(1, 2),
                              " Degree of the pressure finite element space. ");
        }
        prm.leave_subsection();

        prm.enter_subsection("Physical coefficients");
        {
            prm.enter_subsection("Viscosity");
            {
                prm.declare_entry("coeff_nu",
                                  "1.825e-5",
                                  Patterns::Double(0.0),
                                  " Viscosity coefficient (nu) [Pa*s]. ");
                prm.declare_entry("coeff_nu_start",
                                  "2.5e-3",
                                  Patterns::Double(0.0),
                                  " Viscosity coefficient start for ramp-down procedure (nu_0) [Pa*s]. ");
                prm.declare_entry("coeff_nu_ramp_down_times",
                                  "10",
                                  Patterns::Integer(1),
                                  " Amount of iterations for ramp-down procedure. ");
            }
            prm.leave_subsection();

            prm.declare_entry("coeff_rho",
                              "1.204",
                              Patterns::Double(0.0),
                              " Density of the fluid [kg/m^3]. ");
        }
        prm.leave_subsection();

        prm.enter_subsection("Inlet velocity");
        {
            prm.declare_entry("inlet_velocity_start",
                              "0.0, 0.0, 0.0",
                              Patterns::List(Patterns::Double(Patterns::Double::min_double_value), 3, 3),
                              " Inlet velocity components at time 0 [m/s]. ");
            prm.declare_entry("inlet_velocity_end",
                              "0.0, 1.0, 0.0",
                              Patterns::List(Patterns::Double(Patterns::Double::min_double_value), 3, 3),
                              " Inlet velocity components at end time [m/s]. ");
        }
        prm.leave_subsection();

        prm.declare_entry("outlet_pressure",
                          "10.0",
                          Patterns::Double(Patterns::Double::min_double_value),
                          " Outlet pressure [Pa]. ");

        prm.enter_subsection("Solver settings");
        {
            prm.declare_entry("max_solver_iteration_amt",
                              "5000",
                              Patterns::Integer(1),
                              " Maximum amount of iterations for the linear solver. ");
            prm.declare_entry("desired_solver_precision",
                              "1e-6",
                              Patterns::Double(0.0),
                              " Desired precision for the linear solver. ");
            prm.declare_entry("max_newton_iteration_amt",
                              "1000",
                              Patterns::Integer(1),
                              " Maximum amount of iterations for the Newton solver. ");
            prm.declare_entry("desired_newton_precision",
                              "1e-6",
                              Patterns::Double(0.0),
                              " Desired precision for the Newton solver. ");
            /*prm.declare_entry("theta",
                              "0.5",
                              Patterns::Double(0.0, 1.0),
                              " Theta parameter for the time discretization. ");*/
            /*prm.declare_entry("coeff_relax_gamma",
                              "1.0",
                              Patterns::Double(0.0),
                              " Relaxation parameter for Netwon. ");*/
            prm.declare_entry("preconditioner_coeff_alpha",
                              "0.1",
                              Patterns::Double(0.0),
                              " Used for the preconditioner. ");
        }
        prm.leave_subsection();

        prm.enter_subsection("Time step data");
        {
            prm.declare_entry("time_steps_pre_ramp",
                              "0",
                              Patterns::Integer(0),
                              " Amount of time steps in which the velocity is stationary, before the ramp starts. ");
            prm.declare_entry("time_steps_ramp",
                              "20",
                              Patterns::Integer(1),
                              " Amount of time steps needed to go from the starting velocity to the ending velocity. ");
            prm.declare_entry("time_steps_post_ramp",
                              "10",
                              Patterns::Integer(0),
                              " Amount of time steps in which the velocity is stationary, after the ramp has ended. ");
            prm.declare_entry("time_steps_per_second",
                              "10",
                              Patterns::Integer(1),
                              " Amount of time steps per second. ");
        }
        prm.leave_subsection();

        prm.enter_subsection("Boundary conditions");
        {
            prm.declare_entry("surfaces_walls",
                              "1, 5",
                              Patterns::List(Patterns::Integer(0)),
                              " Wall surfaces indices ");
            prm.declare_entry("surfaces_inlets",
                              "2",
                              Patterns::List(Patterns::Integer(0)),
                              " Inlet surfaces indices. ");
            prm.declare_entry("surfaces_outlets",
                              "3",
                              Patterns::List(Patterns::Integer(0)),
                              " Inlet surfaces indices. ");
            prm.declare_entry("surfaces_free_slip",
                              "4",
                              Patterns::List(Patterns::Integer(0)),
                              " Free-slip surfaces indices. ");
        }
        prm.leave_subsection();

        prm.declare_entry("characteristic_length",
                          "0.3",
                          Patterns::Double(0.0),
                          " Characteristic length used to estimate Reynold's number. ");

        update_variables();
    }

    void SimulationSettings::print() const
    {
        prm.print_parameters(std::cout, ParameterHandler::OutputStyle::ShortJSON);
        // prm.print_parameters("params.json");
    }

    void print_simulation_settings(const SimulationSettings &sim_settings)
    {
        sim_settings.print();
    }

    void SimulationSettings::update_variables()
    {
        file_name = prm.get("file_name");
        out_file_name = prm.get("out_file_name");

        prm.enter_subsection("Finite-element settings");
        {
            degree_velocity = prm.get_integer("degree_velocity");
            degree_pressure = prm.get_integer("degree_pressure");
        }
        prm.leave_subsection();

        prm.enter_subsection("Physical coefficients");
        {
            prm.enter_subsection("Viscosity");
            {
                coeff_nu = prm.get_double("coeff_nu");
                coeff_nu_start = prm.get_double("coeff_nu_start");
                coeff_nu_ramp_down_times = prm.get_integer("coeff_nu_ramp_down_times");
            }
            prm.leave_subsection();

            coeff_rho = prm.get_double("coeff_rho");
        }
        prm.leave_subsection();

        prm.enter_subsection("Inlet velocity");
        {
            std::istringstream temp;

            temp.str(prm.get("inlet_velocity_start"));

            for (unsigned int i = 0; i < 3; ++i)
            {
                temp >> inlet_velocity_start[i];
                if (inlet_velocity_start[i] < 1e-50)
                    inlet_velocity_start[i] = 0.0;
            }

            temp.clear();

            temp.str(prm.get("inlet_velocity_end"));

            for (unsigned int i = 0; i < 3; ++i)
            {
                temp >> inlet_velocity_end[i];
                if (inlet_velocity_end[i] < 1e-50)
                    inlet_velocity_end[i] = 0.0;
            }

            temp.clear();
        }
        prm.leave_subsection();

        outlet_pressure = prm.get_double("outlet_pressure");

        prm.enter_subsection("Solver settings");
        {
            max_solver_iteration_amt = prm.get_integer("max_solver_iteration_amt");
            desired_solver_precision = prm.get_double("desired_solver_precision");

            max_newton_iteration_amt = prm.get_integer("max_newton_iteration_amt");
            desired_newton_precision = prm.get_double("desired_newton_precision");

            // theta = prm.get_double("theta");
            preconditioner_coeff_alpha = prm.get_double("preconditioner_coeff_alpha");
            // coeff_relax_gamma = prm.get_double("coeff_relax_gamma");
        }
        prm.leave_subsection();

        prm.enter_subsection("Time step data");
        {
            time_steps_pre_ramp = prm.get_integer("time_steps_pre_ramp");
            time_steps_ramp = prm.get_integer("time_steps_ramp");
            time_steps_post_ramp = prm.get_integer("time_steps_post_ramp");
            time_steps_per_second = prm.get_double("time_steps_per_second");
        }
        prm.leave_subsection();

        prm.enter_subsection("Boundary conditions");
        {
            std::istringstream temp;
            unsigned int temp_uint;

            temp.str(prm.get("surfaces_walls"));

            surfaces_walls.clear();

            while (temp >> temp_uint)
            {
                surfaces_walls.push_back(temp_uint);
            }
            temp.clear();

            temp.str(prm.get("surfaces_inlets"));

            surfaces_inlets.clear();
            while (temp >> temp_uint)
            {
                surfaces_inlets.push_back(temp_uint);
            }

            temp.clear();

            temp.str(prm.get("surfaces_outlets"));

            surfaces_outlets.clear();
            while (temp >> temp_uint)
            {
                surfaces_outlets.push_back(temp_uint);
            }

            temp.clear();

            temp.str(prm.get("surfaces_free_slip"));

            surfaces_free_slip.clear();
            while (temp >> temp_uint)
            {
                surfaces_free_slip.push_back(temp_uint);
            }

            temp.clear();
        }
        prm.leave_subsection();

        characteristic_length = prm.get_double("characteristic_length");

        int mpi_ID;
        int mpi_size;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_ID);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        bool warning_issued = false;
        if (mpi_ID == 0)
        {
            if (surfaces_inlets.size() == 0)
            {
                std::cout << "Warning: no inlet surfaces specified." << std::endl;
                warning_issued = true;
            }
            if (surfaces_outlets.size() == 0)
            {
                std::cout << "Warning: no outlet surfaces specified." << std::endl;
                warning_issued = true;
            }
            if (surfaces_free_slip.size() == 0)
            {
                std::cout << "Warning: no free-slip surfaces specified." << std::endl;
                warning_issued = true;
            }
            if (surfaces_walls.size() == 0)
            {
                std::cout << "Warning: no wall surfaces specified." << std::endl;
                warning_issued = true;
            }

            if (warning_issued)
            {
                std::cout << "Please check whether the input syntax is correct." << std::endl;
            }
        }
    }

    void SimulationSettings::read_data(const std::string &filename)
    {
        std::ifstream file(filename, std::ios::in);
        AssertThrow(file, ExcFileNotOpen(filename));

        prm.parse_input_from_json(file);

        update_variables();

        file.close();
    }
}

/*
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

*/