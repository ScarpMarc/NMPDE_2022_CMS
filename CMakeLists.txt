cmake_minimum_required(VERSION 3.12.0)
project(PteroWindChamber LANGUAGES CXX C)

#execute_process (
#    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/Modules/load_modules.sh
#    COMMAND_ECHO STDOUT
#)

include(common/cmake-common.cmake)

add_executable(Ptero src/main.cpp src/Problem.cpp src/SimulationSettings.cpp)

#add_custom_command(
#  TARGET Ptero PRE_BUILD
#  COMMAND module load gcc-glibc dealii
#  )
#add_custom_command(
#  TARGET Ptero POST_BUILD
#  COMMAND module purge
#  )


add_custom_command(
  TARGET Ptero POST_BUILD
  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/mesh/Ptero_full.msh ${CMAKE_CURRENT_BINARY_DIR}/
  )

deal_ii_setup_target(Ptero)

#execute_process (
#    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/Modules/unload_modules.sh
#    COMMAND_ECHO STDOUT
#)
