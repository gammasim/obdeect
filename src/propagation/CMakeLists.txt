cmake_minimum_required(VERSION 3.20)
project(obdeect_propagation LANGUAGES CXX)

option(OBDEECT_PROPAGATION_BUILD_TESTING "Build the C++ tests" ON)
if(DEFINED SKBUILD_PLATLIB_DIR)
  set(OBDEECT_PROPAGATION_BUILD_PYTHON_API_DEFAULT ON)
else()
  set(OBDEECT_PROPAGATION_BUILD_PYTHON_API_DEFAULT OFF)
endif()
option(OBDEECT_PROPAGATION_BUILD_PYTHON_API "Build the Python extension" ${OBDEECT_PROPAGATION_BUILD_PYTHON_API_DEFAULT})
option(OBDEECT_PROPAGATION_BUILD_TOOLS "Build the obdeect propagation run command-line tool" ON)

# The stream component: the target of the combined build, else the stream sources next to
# this component (builds from the repository), else the copy bundled into the sdist
# (stream-src, staged at release time), else an installed package.
if(NOT TARGET obdeect::stream)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../stream/CMakeLists.txt)
    set(OBDEECT_STREAM_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../stream)
  elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/stream-src/CMakeLists.txt)
    set(OBDEECT_STREAM_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/stream-src)
  endif()
  if(OBDEECT_STREAM_SOURCE_DIR)
    set(OBDEECT_STREAM_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(OBDEECT_STREAM_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
    set(OBDEECT_STREAM_BUILD_CAPI OFF CACHE BOOL "" FORCE)
    set(OBDEECT_STREAM_BUILD_PYTHON_API OFF CACHE BOOL "" FORCE)
    set(OBDEECT_STREAM_INSTALL OFF CACHE BOOL "" FORCE)
    add_subdirectory(${OBDEECT_STREAM_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/stream)
    message(STATUS "obdeect propagation: using the stream sources in ${OBDEECT_STREAM_SOURCE_DIR}")
  else()
    find_package(Python 3.14 COMPONENTS Interpreter)
    if(Python_Interpreter_FOUND AND NOT obdeect_stream_DIR)
      execute_process(COMMAND "${Python_EXECUTABLE}" -m obdeect.stream --cmake-dir
                      OUTPUT_VARIABLE obdeect_stream_DIR OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    endif()
    find_package(obdeect_stream CONFIG REQUIRED)
    message(STATUS "obdeect propagation: using the installed stream package ${obdeect_stream_VERSION} (${obdeect_stream_DIR})")
  endif()
endif()

# Header-only, standard-library-only core.
add_library(obdeect_propagation_core INTERFACE)
target_include_directories(obdeect_propagation_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include)
target_compile_features(obdeect_propagation_core INTERFACE cxx_std_20)
target_link_libraries(obdeect_propagation_core INTERFACE obdeect::stream)
# Identical floating-point results for any thread count and chunk size on one
# platform: no implicit fused multiply-add contraction.
target_compile_options(obdeect_propagation_core INTERFACE -ffp-contract=off)

function(obdeect_propagation_enable_warnings target)
  target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
endfunction()

if(OBDEECT_PROPAGATION_BUILD_PYTHON_API)
  find_package(Python 3.14 COMPONENTS Interpreter Development.Module REQUIRED)
  execute_process(COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
                  OUTPUT_VARIABLE nanobind_DIR OUTPUT_STRIP_TRAILING_WHITESPACE
                  COMMAND_ERROR_IS_FATAL ANY)
  find_package(nanobind CONFIG REQUIRED)
  nanobind_add_module(_core NB_STATIC cpp/bindings/core.cpp)
  target_link_libraries(_core PRIVATE obdeect_propagation_core)
  obdeect_propagation_enable_warnings(_core)
  install(TARGETS _core LIBRARY DESTINATION obdeect/propagation)
endif()

if(OBDEECT_PROPAGATION_BUILD_TOOLS)
  add_executable(obdeect-propagation-run cpp/tools/propagate.cpp)
  target_link_libraries(obdeect-propagation-run PRIVATE obdeect_propagation_core)
  obdeect_propagation_enable_warnings(obdeect-propagation-run)
  if(OBDEECT_PROPAGATION_BUILD_TESTING)
    add_test(NAME run_version COMMAND obdeect-propagation-run --version)
    set_tests_properties(run_version PROPERTIES PASS_REGULAR_EXPRESSION
      "^obdeect propagation run [0-9]+\\.[0-9]+\\.[0-9]+ \\(photon-stream format [0-9]+\\.[0-9]+\\)")
  endif()
  # Wheels: inside the package, run through `obdeect propagation run`.
  if(DEFINED SKBUILD_PLATLIB_DIR)
    install(TARGETS obdeect-propagation-run RUNTIME DESTINATION obdeect/propagation/bin)
  endif()
endif()

if(OBDEECT_PROPAGATION_BUILD_TESTING)
  enable_testing()
  foreach(test_name version transmission thinning pipeline profile refraction timing
                    propagator scattering stream_contract physics)
    add_executable(obdeect_propagation_test_${test_name} cpp/tests/test_${test_name}.cpp)
    target_link_libraries(obdeect_propagation_test_${test_name} PRIVATE obdeect_propagation_core)
    # check.hpp lives with the stream tests; tests are built from the repository only, never from an sdist.
    target_include_directories(obdeect_propagation_test_${test_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../stream/cpp/tests)
    # Reference data: cpp/tests/data (see test_transmission.cpp).
    target_compile_definitions(obdeect_propagation_test_${test_name}
                               PRIVATE OBDEECT_PROPAGATION_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/cpp/tests/data")
    obdeect_propagation_enable_warnings(obdeect_propagation_test_${test_name})
    add_test(NAME ${test_name} COMMAND obdeect_propagation_test_${test_name})
  endforeach()
endif()
