project(dagmc_model_reporter)

cmake_minimum_required(VERSION 3.10)

find_package(DAGMC REQUIRED)

add_subdirectory(vendor/argparse)

add_executable(dagmc-reporter ${CMAKE_CURRENT_SOURCE_DIR}/reporter.cpp)
target_link_libraries(dagmc-reporter dagmc-shared argparse)
