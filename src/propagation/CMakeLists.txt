cmake_minimum_required(VERSION 3.20)
project(obdeect_propagation LANGUAGES CXX)

# Shared CMake helpers: the repository's cmake/ directory, or the copy staged next to the
# sdist (shared-cmake/). A parent that embeds this component has set the directory already.
if(NOT DEFINED OBDEECT_SHARED_CMAKE_DIR)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/ObdeectCommon.cmake)
    set(OBDEECT_SHARED_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake)
  else()
    set(OBDEECT_SHARED_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/shared-cmake)
  endif()
endif()
include(${OBDEECT_SHARED_CMAKE_DIR}/ObdeectCommon.cmake)

option(OBDEECT_PROPAGATION_BUILD_TESTING "Build the C++ tests" ON)
if(DEFINED SKBUILD_PLATLIB_DIR)
  set(OBDEECT_PROPAGATION_BUILD_PYTHON_API_DEFAULT ON)
else()
  set(OBDEECT_PROPAGATION_BUILD_PYTHON_API_DEFAULT OFF)
endif()
option(OBDEECT_PROPAGATION_BUILD_PYTHON_API "Build the Python extension" ${OBDEECT_PROPAGATION_BUILD_PYTHON_API_DEFAULT})
option(OBDEECT_PROPAGATION_BUILD_TOOLS "Build the obdeect propagation run command-line tool" ON)

# The stream component (see obdeect_locate_stream in cmake/ObdeectCommon.cmake).
obdeect_locate_stream(propagation)

# Header-only, standard-library-only core.
add_library(obdeect_propagation_core INTERFACE)
target_include_directories(obdeect_propagation_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include)
find_package(Threads REQUIRED)
target_link_libraries(obdeect_propagation_core INTERFACE obdeect::stream Threads::Threads)
obdeect_set_core_options(obdeect_propagation_core)

if(OBDEECT_PROPAGATION_BUILD_PYTHON_API)
  obdeect_add_python_module(propagation obdeect_propagation_core)
endif()

if(OBDEECT_PROPAGATION_BUILD_TOOLS)
  add_executable(obdeect-propagation-run cpp/tools/propagate.cpp)
  target_link_libraries(obdeect-propagation-run PRIVATE obdeect_propagation_core)
  obdeect_enable_warnings(obdeect-propagation-run)
  if(OBDEECT_PROPAGATION_BUILD_TESTING)
    add_test(NAME run_version COMMAND obdeect-propagation-run --version)
    set_tests_properties(run_version PROPERTIES PASS_REGULAR_EXPRESSION
      "^obdeect propagation run [0-9]+\\.[0-9]+\\.[0-9]+ \\(photon-stream format [0-9]+\\.[0-9]+\\)")
  endif()
  # Wheels: inside the package, run through `obdeect propagation run`.
  if(DEFINED SKBUILD_PLATLIB_DIR)
    install(TARGETS obdeect-propagation-run RUNTIME DESTINATION obdeect/propagation/bin)
  endif()
  add_executable(obdeect-propagation-to-sky cpp/tools/to_sky.cpp)
  target_link_libraries(obdeect-propagation-to-sky PRIVATE obdeect_propagation_core)
  obdeect_enable_warnings(obdeect-propagation-to-sky)
  if(OBDEECT_PROPAGATION_BUILD_TESTING)
    add_test(NAME to_sky_version COMMAND obdeect-propagation-to-sky --version)
    set_tests_properties(to_sky_version PROPERTIES PASS_REGULAR_EXPRESSION
      "^obdeect propagation to-sky [0-9]+\\.[0-9]+\\.[0-9]+ \\(photon-stream format [0-9]+\\.[0-9]+\\)")
  endif()
  if(DEFINED SKBUILD_PLATLIB_DIR)
    install(TARGETS obdeect-propagation-to-sky RUNTIME DESTINATION obdeect/propagation/bin)
  endif()
endif()

if(OBDEECT_PROPAGATION_BUILD_TESTING)
  enable_testing()
  # Tests that take more than about 5 s (labelled slow, see cmake/ObdeectTests.cmake).
  set(OBDEECT_PROPAGATION_SLOW_TESTS curvature_modes dispersion scattering sphere
      curved_refraction propagator)
  foreach(test_name version transmission thinning pipeline profile refraction timing
                    propagator scattering stream_contract physics curvature curvature_modes telescopes
                    curved_refraction dispersion horizon to_sky upward sphere chord)
    add_executable(obdeect_propagation_test_${test_name} cpp/tests/test_${test_name}.cpp)
    target_link_libraries(obdeect_propagation_test_${test_name} PRIVATE obdeect_propagation_core)
    # check.hpp lives with the stream tests; tests are built from the repository only, never from an sdist.
    target_include_directories(obdeect_propagation_test_${test_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../stream/cpp/tests)
    # Reference data: cpp/tests/data (see test_transmission.cpp).
    target_compile_definitions(obdeect_propagation_test_${test_name}
                               PRIVATE OBDEECT_PROPAGATION_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/cpp/tests/data")
    obdeect_enable_warnings(obdeect_propagation_test_${test_name})
    set(slow_arg)
    if(test_name IN_LIST OBDEECT_PROPAGATION_SLOW_TESTS)
      set(slow_arg SLOW)
    endif()
    obdeect_add_test(${test_name} propagation SOURCE test_${test_name}.cpp ${slow_arg}
                     COMMAND obdeect_propagation_test_${test_name})
  endforeach()
  obdeect_check_test_sources(${CMAKE_CURRENT_SOURCE_DIR}/cpp/tests)
endif()
