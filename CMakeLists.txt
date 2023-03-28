cmake_minimum_required(VERSION 3.12.0)
project(PteroWindChamber LANGUAGES CXX C)

include(common/cmake-common.cmake)

add_executable(Ptero src/main.cpp src/Problem.cpp src/SimulationSettings.cpp)
deal_ii_setup_target(Ptero)
