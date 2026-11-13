cmake_minimum_required(VERSION 3.20)
project(dlufd VERSION 0.1.0 LANGUAGES CXX)

option(DLUFD_BUILD_TESTING "Build the C++ tests" ON)
if(DEFINED SKBUILD_PLATLIB_DIR)
  set(DLUFD_BUILD_PYTHON_API_DEFAULT ON)
else()
  set(DLUFD_BUILD_PYTHON_API_DEFAULT OFF)
endif()
option(DLUFD_BUILD_PYTHON_API "Build the Python extension" ${DLUFD_BUILD_PYTHON_API_DEFAULT})
option(DLUFD_BUILD_TOOLS "Build the dlufd-propagate command-line tool" ON)

# The stream component: the target of the combined build, else the stream sources next to
# this component (builds from the repository), else an installed package (sdist builds).
if(NOT TARGET gstrickd::gstrickd_core)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../stream/CMakeLists.txt)
    set(GSTRICKD_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(GSTRICKD_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
    set(GSTRICKD_BUILD_CAPI OFF CACHE BOOL "" FORCE)
    set(GSTRICKD_BUILD_PYTHON_API OFF CACHE BOOL "" FORCE)
    set(GSTRICKD_INSTALL OFF CACHE BOOL "" FORCE)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../stream ${CMAKE_CURRENT_BINARY_DIR}/stream)
    message(STATUS "dlufd: using the stream sources in ../stream")
  else()
    find_package(Python 3.14 COMPONENTS Interpreter)
    if(Python_Interpreter_FOUND AND NOT gstrickd_DIR)
      execute_process(COMMAND "${Python_EXECUTABLE}" -m gstrickd --cmake-dir
                      OUTPUT_VARIABLE gstrickd_DIR OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    endif()
    find_package(gstrickd CONFIG REQUIRED)
    message(STATUS "dlufd: using installed gstrickd ${gstrickd_VERSION} (${gstrickd_DIR})")
  endif()
endif()

# Header-only, standard-library-only core.
add_library(dlufd_core INTERFACE)
target_include_directories(dlufd_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include)
target_compile_features(dlufd_core INTERFACE cxx_std_20)
target_link_libraries(dlufd_core INTERFACE gstrickd::gstrickd_core)
# Identical floating-point results for any thread count and chunk size on one
# platform: no implicit fused multiply-add contraction.
target_compile_options(dlufd_core INTERFACE -ffp-contract=off)

function(dlufd_enable_warnings target)
  target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
endfunction()

if(DLUFD_BUILD_PYTHON_API)
  find_package(Python 3.14 COMPONENTS Interpreter Development.Module REQUIRED)
  execute_process(COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
                  OUTPUT_VARIABLE nanobind_DIR OUTPUT_STRIP_TRAILING_WHITESPACE
                  COMMAND_ERROR_IS_FATAL ANY)
  find_package(nanobind CONFIG REQUIRED)
  nanobind_add_module(_core NB_STATIC cpp/bindings/core.cpp)
  target_link_libraries(_core PRIVATE dlufd_core)
  dlufd_enable_warnings(_core)
  install(TARGETS _core LIBRARY DESTINATION dlufd)
endif()

if(DLUFD_BUILD_TOOLS)
  add_executable(dlufd-propagate cpp/tools/dlufd_propagate.cpp)
  target_link_libraries(dlufd-propagate PRIVATE dlufd_core)
  dlufd_enable_warnings(dlufd-propagate)
  if(DEFINED SKBUILD_SCRIPTS_DIR)
    install(TARGETS dlufd-propagate RUNTIME DESTINATION ${SKBUILD_SCRIPTS_DIR})
  endif()
endif()

if(DLUFD_BUILD_TESTING)
  enable_testing()
  foreach(test_name version transmission thinning pipeline profile refraction timing
                    propagator scattering gstrickd_contract physics)
    add_executable(dlufd_test_${test_name} cpp/tests/test_${test_name}.cpp)
    target_link_libraries(dlufd_test_${test_name} PRIVATE dlufd_core)
    # Reference data: cpp/tests/data (see test_transmission.cpp).
    target_compile_definitions(dlufd_test_${test_name}
                               PRIVATE DLUFD_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/cpp/tests/data")
    dlufd_enable_warnings(dlufd_test_${test_name})
    add_test(NAME ${test_name} COMMAND dlufd_test_${test_name})
  endforeach()
endif()
