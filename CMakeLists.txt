cmake_minimum_required(VERSION 3.12.0)

# Locate deal.II and initialize its variables.
find_package(deal.II 9.3.1 REQUIRED
  HINTS ${DEAL_II_DIR} $ENV{DEAL_II_DIR} $ENV{mkDealiiPrefix})
deal_ii_initialize_cached_variables()

project(PteroWindChamber LANGUAGES CXX C)

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


#add_custom_command(
#  TARGET Ptero POST_BUILD
#  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/mesh/* ${CMAKE_CURRENT_SOURCE_DIR}/input/* ${CMAKE_CURRENT_BINARY_DIR}/
#  )

deal_ii_setup_target(Ptero)

add_custom_target("deploy")
add_custom_command(TARGET "deploy" POST_BUILD 
COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/mesh/* ${CMAKE_CURRENT_SOURCE_DIR}/input/* ${CMAKE_CURRENT_BINARY_DIR}/
  )

#execute_process (
#    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/Modules/unload_modules.sh
#    COMMAND_ECHO STDOUT
#)
