cmake_minimum_required(VERSION 3.20)
project(solver VERSION 0.1.0 LANGUAGES CXX)

option(SOLVER_BUILD_TESTS "Build solver tests" ON)
option(SOLVER_BUILD_EXAMPLES "Build solver examples" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)
find_package(Eigen3 3.4 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
  set(_solver_build_testing_was_defined FALSE)
  if(DEFINED BUILD_TESTING)
    set(_solver_build_testing_was_defined TRUE)
    set(_solver_saved_build_testing "${BUILD_TESTING}")
  endif()

  set(BUILD_TESTING OFF CACHE BOOL "Disable dependency tests while fetching Eigen" FORCE)
  FetchContent_Declare(
    eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
  FetchContent_MakeAvailable(eigen)

  if(_solver_build_testing_was_defined)
    set(BUILD_TESTING "${_solver_saved_build_testing}" CACHE BOOL "Build tests" FORCE)
  else()
    unset(BUILD_TESTING CACHE)
  endif()
endif()

add_library(solver_core
  src/core/DofManager.cpp
  src/core/Model.cpp
  src/elements/Beam3D.cpp
  src/solver/LinearStaticSolver.cpp
)
add_library(solver::core ALIAS solver_core)

target_include_directories(solver_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(solver_core PUBLIC Eigen3::Eigen)

target_compile_options(solver_core PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra;-Wpedantic;-Wconversion>
  $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

if(SOLVER_BUILD_EXAMPLES)
  add_executable(cantilever_example examples/cantilever.cpp)
  target_link_libraries(cantilever_example PRIVATE solver::core)
endif()

if(SOLVER_BUILD_TESTS)
  enable_testing()
  add_executable(test_beam3d tests/test_beam3d.cpp)
  target_link_libraries(test_beam3d PRIVATE solver::core)
  add_test(NAME beam3d COMMAND test_beam3d)
endif()
