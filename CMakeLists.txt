cmake_minimum_required(VERSION 3.20)
project(solver VERSION 0.2.0 LANGUAGES CXX)

option(SOLVER_BUILD_TESTS "Build solver tests" ON)
option(SOLVER_BUILD_EXAMPLES "Build solver examples" ON)
option(SOLVER_BUILD_BENCHMARKS "Build solver performance benchmarks" ON)

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

find_package(Spectra 1.2 QUIET)
if(NOT TARGET Spectra::Spectra)
  FetchContent_Declare(
    spectra
    URL https://github.com/yixuan/spectra/archive/refs/tags/v1.2.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
  FetchContent_MakeAvailable(spectra)
  if(TARGET Spectra AND NOT TARGET Spectra::Spectra)
    add_library(Spectra::Spectra ALIAS Spectra)
  endif()
endif()

add_library(solver_core
  src/core/DofManager.cpp
  src/core/Model.cpp
  src/dynamics/RayleighDamping.cpp
  src/elements/Beam3D.cpp
  src/io/BeamForceCsvWriter.cpp
  src/io/BeamForceDiagramSvgWriter.cpp
  src/loads/BeamLinearLoad3D.cpp
  src/loads/BeamPointLoad3D.cpp
  src/loads/BeamPartialLinearLoad3D.cpp
  src/loads/BeamUniformLoad3D.cpp
  src/loads/TimeDependentElementLoad.cpp
  src/recorder/ElementRecorder.cpp
  src/recorder/EnvelopeRecorder.cpp
  src/recorder/NodeRecorder.cpp
  src/recorder/SectionRecorder.cpp
  src/response/BeamResponseSampler.cpp
  src/response/BeamSectionStressRecovery.cpp
  src/solver/LinearStaticSolver.cpp
  src/solver/ModalSolver.cpp
  src/solver/NewmarkBetaSolver.cpp
  src/solver/SparseDofReducer.cpp
)
add_library(solver::core ALIAS solver_core)

target_include_directories(solver_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(solver_core
  PUBLIC Eigen3::Eigen
  PRIVATE Spectra::Spectra
)

target_compile_options(solver_core PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra;-Wpedantic;-Wconversion>
  $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

if(SOLVER_BUILD_EXAMPLES)
  add_executable(cantilever_example examples/cantilever.cpp)
  target_link_libraries(cantilever_example PRIVATE solver::core)

  add_executable(dynamic_cantilever_example examples/dynamic_cantilever.cpp)
  target_link_libraries(dynamic_cantilever_example PRIVATE solver::core)

  add_executable(beam_force_diagrams_example examples/beam_force_diagrams.cpp)
  target_link_libraries(beam_force_diagrams_example PRIVATE solver::core)

  add_executable(beam_stress_recovery_example examples/beam_stress_recovery.cpp)
  target_link_libraries(beam_stress_recovery_example PRIVATE solver::core)

  add_executable(result_envelope_example examples/result_envelope.cpp)
  target_link_libraries(result_envelope_example PRIVATE solver::core)
endif()

if(SOLVER_BUILD_BENCHMARKS)
  add_executable(solver_benchmark benchmarks/solver_benchmark.cpp)
  target_link_libraries(solver_benchmark PRIVATE solver::core)
endif()

if(SOLVER_BUILD_TESTS)
  enable_testing()

  add_executable(test_beam3d tests/test_beam3d.cpp)
  target_link_libraries(test_beam3d PRIVATE solver::core)
  add_test(NAME beam3d COMMAND test_beam3d)

  add_executable(test_dynamics tests/test_dynamics.cpp)
  target_link_libraries(test_dynamics PRIVATE solver::core)
  add_test(NAME dynamics COMMAND test_dynamics)

  add_executable(test_element_response tests/test_element_response.cpp)
  target_link_libraries(test_element_response PRIVATE solver::core)
  add_test(NAME element_response COMMAND test_element_response)

  add_executable(test_beam_loads tests/test_beam_loads.cpp)
  target_link_libraries(test_beam_loads PRIVATE solver::core)
  add_test(NAME beam_loads COMMAND test_beam_loads)

  add_executable(test_partial_beam_loads tests/test_partial_beam_loads.cpp)
  target_link_libraries(test_partial_beam_loads PRIVATE solver::core)
  add_test(NAME partial_beam_loads COMMAND test_partial_beam_loads)

  add_executable(test_time_element_loads tests/test_time_element_loads.cpp)
  target_link_libraries(test_time_element_loads PRIVATE solver::core)
  add_test(NAME time_element_loads COMMAND test_time_element_loads)

  add_executable(test_beam_section_forces tests/test_beam_section_forces.cpp)
  target_link_libraries(test_beam_section_forces PRIVATE solver::core)
  add_test(NAME beam_section_forces COMMAND test_beam_section_forces)

  add_executable(test_beam_force_output tests/test_beam_force_output.cpp)
  target_link_libraries(test_beam_force_output PRIVATE solver::core)
  add_test(NAME beam_force_output COMMAND test_beam_force_output)

  add_executable(test_section_stress tests/test_section_stress.cpp)
  target_link_libraries(test_section_stress PRIVATE solver::core)
  add_test(NAME section_stress COMMAND test_section_stress)

  add_executable(test_result_recorders tests/test_result_recorders.cpp)
  target_link_libraries(test_result_recorders PRIVATE solver::core)
  add_test(NAME result_recorders COMMAND test_result_recorders)

  add_executable(test_sparse_assembly tests/test_sparse_assembly.cpp)
  target_link_libraries(test_sparse_assembly PRIVATE solver::core)
  add_test(NAME sparse_assembly COMMAND test_sparse_assembly)

  add_executable(test_sparse_solvers tests/test_sparse_solvers.cpp)
  target_link_libraries(test_sparse_solvers PRIVATE solver::core)
  add_test(NAME sparse_solvers COMMAND test_sparse_solvers)

  add_executable(test_v1_verification tests/test_v1_verification.cpp)
  target_link_libraries(test_v1_verification PRIVATE solver::core)
  add_test(NAME v1_verification COMMAND test_v1_verification)

  add_executable(test_sparse_modal tests/test_sparse_modal.cpp)
  target_link_libraries(test_sparse_modal PRIVATE solver::core)
  add_test(NAME sparse_modal COMMAND test_sparse_modal)
endif()
