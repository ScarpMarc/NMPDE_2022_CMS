#include "SimulationSettings.hpp"
#include <iostream>
#include <fstream>

using namespace dealii;

namespace ns_sim_settings
{
    SimulationSettings::SimulationSettings()
        : SimulationSettings(
              "Ptero_full.msh", // file_name("Ptero_full.msh"),
              2,                //        degree_velocity(2),
              1,                //        degree_pressure(1),
              1.825e-5,         //        coeff_nu(1.825e-5),
              2.5e-3,           //        coeff_nu_start(2.5e-3),
              10,               //        coeff_nu_ramp_down_times(10),
              1.204,            //        coeff_rho(1.204),
              {0.0, 10.0, 0.0}, //        inlet_velocity({0.0, 10.0, 0.0}),
              {0.0, .1, 0.0},   //        inlet_velocity_increment_per_step({0.0, .1, 0.0}),
              10.0,             //        outlet_pressure(10.0),
              5000,             //        max_solver_iteration_amt(5000),
              1e-6,             //        desired_solver_precision(1e-6),
              1000,             //        max_newton_iteration_amt(1000),
              1e-6,             //        desired_newton_precision(1e-6),
              0.5,              //        theta(0.5),
              100,              //        total_time_steps(100),
              10                //        time_steps_per_second(10)
          )
    {
    }

    SimulationSettings::SimulationSettings(const std::string file_name,
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
                       const unsigned long time_steps_per_second) : file_name(file_name),
                                                                    degree_velocity(degree_velocity),
                                                                    degree_pressure(degree_pressure),
                                                                    coeff_nu(coeff_nu),
                                                                    coeff_nu_start(coeff_nu_start),
                                                                    coeff_nu_ramp_down_times(coeff_nu_ramp_down_times),
                                                                    coeff_rho(coeff_rho),
                                                                    inlet_velocity(inlet_velocity),
                                                                    inlet_velocity_increment_per_step(inlet_velocity_increment_per_step),
                                                                    outlet_pressure(outlet_pressure),
                                                                    max_solver_iteration_amt(max_solver_iteration_amt),
                                                                    desired_solver_precision(desired_solver_precision),
                                                                    max_newton_iteration_amt(max_newton_iteration_amt),
                                                                    desired_newton_precision(desired_newton_precision),
                                                                    theta(theta),
                                                                    total_time_steps(total_time_steps),
                                                                    time_steps_per_second(time_steps_per_second),
                                                                    ramp_down_step_size(0.0),
                                                                    advised_ramp_down_step_size(0.0)
    {
        prm.declare_entry("file_name",
                          "Ptero_full.msh",
                          Patterns::FileName(),
                          " Mesh file name.");

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
            prm.declare_entry("inlet_velocity",
                              "0.0, 10.0, 0.0",
                              Patterns::List(Patterns::Double(0.0)),
                              " Inlet velocity components [m/s]. ");
            prm.declare_entry("inlet_velocity_increment_per_step",
                              "0.0, 0.1, 0.0",
                              Patterns::List(Patterns::Double(0.0)),
                              " Inlet velocity increment per time step [m/s]. ");
        }
        prm.leave_subsection();

        prm.declare_entry("outlet_pressure",
                          "10.0",
                          Patterns::Double(0.0),
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
            prm.declare_entry("theta",
                              "0.5",
                              Patterns::Double(0.0, 1.0),
                              " Theta parameter for the time discretization. ");
        }

        prm.enter_subsection("Time step data");
        {
            prm.declare_entry("total_time_steps",
                              "100",
                              Patterns::Integer(1),
                              " Total amount of time steps. ");
            prm.declare_entry("time_steps_per_second",
                              "10",
                              Patterns::Integer(1),
                              " Amount of time steps per second. ");
        }
        prm.leave_subsection();

        prm.enter_subsection("Space discretization");
        {
            prm.declare_entry("n_of_refines",
                              "0",
                              Patterns::Integer(0, 15),
                              " The number of global refines we do on the mesh. ");
            prm.declare_entry("pressure_fe_degree",
                              "1",
                              Patterns::Integer(1, 5),
                              " The polynomial degree for the pressure space. ");
        }
        prm.leave_subsection();

        prm.enter_subsection("Data solve velocity");
        {
            prm.declare_entry(
                "max_iterations",
                "1000",
                Patterns::Integer(1, 1000),
                " The maximal number of iterations GMRES must make. ");
            prm.declare_entry("eps",
                              "1e-12",
                              Patterns::Double(0.),
                              " The stopping criterion. ");
            prm.declare_entry("Krylov_size",
                              "30",
                              Patterns::Integer(1),
                              " The size of the Krylov subspace to be used. ");
            prm.declare_entry("off_diagonals",
                              "60",
                              Patterns::Integer(0),
                              " The number of off-diagonal elements ILU must "
                              "compute. ");
            prm.declare_entry("diag_strength",
                              "0.01",
                              Patterns::Double(0.),
                              " Diagonal strengthening coefficient. ");
            prm.declare_entry("update_prec",
                              "15",
                              Patterns::Integer(1),
                              " This number indicates how often we need to "
                              "update the preconditioner");
        }
        prm.leave_subsection();

        prm.declare_entry("verbose",
                          "true",
                          Patterns::Bool(),
                          " This indicates whether the output of the solution "
                          "process should be verbose. ");

        prm.declare_entry("output_interval",
                          "1",
                          Patterns::Integer(1),
                          " This indicates between how many time steps we print "
                          "the solution. ");
    }

    void SimulationSettings::print() const
    {
        prm.print_parameters(std::cout, ParameterHandler::OutputStyle::PRM);
    }

    void print_simulation_settings(const SimulationSettings &sim_settings)
    {
        // Output settings
        /*std::cout << "File name: " << sim_settings.file_name << std::endl; //"half_ptero_reducedx2.msh",//std::string file_name;
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

        std::cout << "Max. Newton solver iterations: " << sim_settings.max_newton_iteration_amt << std::endl;   // 1000,//unsigned int max_solver_iteration_amt;
        std::cout << "Desired Newton solver precision: " << sim_settings.desired_newton_precision << std::endl; // 1e-6,//double desired_solver_precision;

        std::cout << "Theta for Theta Method: " << sim_settings.theta << std::endl; // 0.5//double theta;
        std::cout << std::endl;
        std::cout << "Total time steps: " << sim_settings.total_time_steps << std::endl;           // 100,// unsigned long total_time_steps;
        std::cout << "Time steps per second: " << sim_settings.time_steps_per_second << std::endl; // 10// double seconds_per_time_steps;
        std::cout << std::endl;*/

        sim_settings.print();
    }

    void SimulationSettings::read_data(const std::string &filename)
    {
        std::ifstream file(filename, std::ios::in);
        AssertThrow(file, ExcFileNotOpen(filename));

        prm.parse_input(file);

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
                coeff_nu = prm.get_double("nu");
                coeff_nu_start = prm.get_double("nu_start");
                coeff_nu_ramp_down_times = prm.get_integer("nu_ramp_down_times");
            }
            prm.leave_subsection();

            coeff_rho = prm.get_double("rho");
        }
        prm.leave_subsection();

        prm.enter_subsection("Inlet velocity");
        {
            prm.enter_subsection("Components");
            {
                inlet_velocity[0] = prm.get_double("x");
                inlet_velocity[1] = prm.get_double("y");
                inlet_velocity[2] = prm.get_double("z");
            }
            prm.leave_subsection();

            prm.enter_subsection("Increment per time step");
            {
                inlet_velocity_increment_per_step[0] = prm.get_double("x");
                inlet_velocity_increment_per_step[1] = prm.get_double("y");
                inlet_velocity_increment_per_step[2] = prm.get_double("z");
            }
            prm.leave_subsection();
        }
        prm.leave_subsection();

        prm.enter_subsection("Outlet pressure");
        {
            outlet_pressure = prm.get_double("p");
        }

        prm.enter_subsection("Solver settings");
        {
            max_solver_iteration_amt = prm.get_integer("max_solver_iteration_amt");
            desired_solver_precision = prm.get_double("desired_solver_precision");

            max_newton_iteration_amt = prm.get_integer("max_newton_iteration_amt");
            desired_newton_precision = prm.get_double("desired_newton_precision");

            theta = prm.get_double("theta");
        }
        prm.leave_subsection();

        prm.enter_subsection("Time settings");
        {
            total_time_steps = prm.get_integer("total_time_steps");
            time_steps_per_second = prm.get_double("time_steps_per_second");
        }
        prm.leave_subsection();

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