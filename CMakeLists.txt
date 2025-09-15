cmake_minimum_required(VERSION 3.16)
project(fmx LANGUAGES CXX)

option(FMX_ENABLE_OPENMP "Enable OpenMP for parallelization" OFF)
option(FMX_SENTMAN_CLOSED_FORM "Use closed-form Sentman expressions" ON)
option(FMX_ENABLE_EMBREE "Enable Embree occlusion backend" OFF)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_library(fmx_core INTERFACE)
target_include_directories(fmx_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(fmx_gsi
  gsi/Sentman.cpp
  gsi/Sentman.hpp
  gsi/CLL.cpp
  gsi/CLL.hpp
  gsi/KernelSet.cpp
  gsi/KernelSet.hpp
  gsi/CLLRuntime.cpp
  gsi/CLLRuntime.hpp
)
target_link_libraries(fmx_gsi PUBLIC fmx_core)
target_include_directories(fmx_gsi PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(FMX_SENTMAN_CLOSED_FORM)
  target_compile_definitions(fmx_gsi PUBLIC FMX_USE_SENTMAN_CLOSED_FORM=1)
endif()

add_library(fmx_geom
  geom/Mesh.cpp
  geom/Mesh.hpp
  geom/Occluder.hpp
  geom/BVH.cpp
  geom/BVH.hpp
)
target_link_libraries(fmx_geom PUBLIC fmx_core)
target_include_directories(fmx_geom PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(fmx_solver
  solver/PanelSolver.cpp
  solver/PanelSolver.hpp
  solver/RegimeAdapter.cpp
  solver/RegimeAdapter.hpp
)
target_link_libraries(fmx_solver PUBLIC fmx_core fmx_gsi fmx_geom)
target_include_directories(fmx_solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

option(FMX_WITH_NRLMSIS2 "Link against NRLMSIS2 if available" OFF)
option(FMX_WITH_HWM14 "Link against HWM14 if available" OFF)
option(FMX_USE_LOCAL_MODELS "Build NRLMSIS2.1 and HWM14 from atm/models with Fortran" OFF)

add_library(fmx_atm
  atm/Atmosphere.hpp
  atm/StubAtmosphere.cpp
  atm/NRLMSIS2.hpp
  atm/NRLMSIS2.cpp
  atm/HWM14.hpp
  atm/HWM14.cpp
  atm/Combined.hpp
)
target_link_libraries(fmx_atm PUBLIC fmx_core)
target_include_directories(fmx_atm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if(FMX_USE_LOCAL_MODELS)
  enable_language(Fortran)
  # NRLMSIS2.1 sources
  file(GLOB NRLMSIS_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/atm/models/NRLMSIS2.1/*.F90
  )
  add_library(nrlmsis2_local STATIC
    ${NRLMSIS_SRC}
    atm/wrappers/nrlmsis_cwrap.f90
  )
  # HWM14 sources
  add_library(hwm14_local STATIC
    atm/models/hwm14/src/hwm14.f90
    atm/models/hwm14/src/hwm_interface.f90
    atm/wrappers/hwm14_cwrap.f90
  )
  target_link_libraries(fmx_atm PUBLIC nrlmsis2_local hwm14_local)
  target_compile_definitions(fmx_atm PUBLIC FMX_WITH_NRLMSIS2_LOCAL=1 FMX_WITH_HWM14_LOCAL=1)
endif()

add_executable(fmx_cli
  cli/fmx.cpp
)
target_link_libraries(fmx_cli PRIVATE fmx_core fmx_gsi fmx_geom fmx_solver fmx_atm)

if(FMX_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(fmx_solver PUBLIC OpenMP::OpenMP_CXX)
  target_link_libraries(fmx_cli PRIVATE OpenMP::OpenMP_CXX)
  target_compile_definitions(fmx_solver PUBLIC FMX_USE_OPENMP=1)
endif()

enable_testing()
add_executable(test_sentman tests/test_sentman.cpp)
target_link_libraries(test_sentman PRIVATE fmx_core fmx_gsi)
add_test(NAME sentman_basic COMMAND test_sentman)
add_executable(test_occlusion tests/test_occlusion.cpp)
target_link_libraries(test_occlusion PRIVATE fmx_core fmx_geom fmx_solver fmx_atm)
add_test(NAME occlusion_basic COMMAND test_occlusion)
add_executable(test_validation tests/test_validation.cpp)
target_link_libraries(test_validation PRIVATE fmx_core fmx_geom fmx_solver fmx_atm)
add_test(NAME validation_trends COMMAND test_validation)
add_executable(test_torque tests/test_torque.cpp)
target_link_libraries(test_torque PRIVATE fmx_core fmx_geom fmx_solver fmx_atm)
add_test(NAME torque_plate_offset COMMAND test_torque)
add_executable(test_cube tests/test_cube.cpp)
target_link_libraries(test_cube PRIVATE fmx_core fmx_geom fmx_solver fmx_atm)
add_test(NAME cube_symmetry_drag COMMAND test_cube)
add_executable(test_cr313 tests/test_cr313.cpp)
target_link_libraries(test_cr313 PRIVATE fmx_core fmx_geom fmx_solver fmx_atm)
add_test(NAME cr313_reference COMMAND test_cr313)
add_executable(test_cll_runtime tests/test_cll_runtime.cpp)
target_link_libraries(test_cll_runtime PRIVATE fmx_core fmx_gsi)
add_test(NAME cll_runtime_basic COMMAND test_cll_runtime)
add_executable(test_cd_cube_dsmc tests/test_cd_cube_dsmc.cpp)
target_link_libraries(test_cd_cube_dsmc PRIVATE fmx_core fmx_geom fmx_solver fmx_atm)
add_test(NAME cd_cube_dsmc_check COMMAND test_cd_cube_dsmc)

add_executable(gen_gsi_table tools/gen_gsi_table.cpp)
target_link_libraries(gen_gsi_table PRIVATE fmx_core fmx_gsi)
