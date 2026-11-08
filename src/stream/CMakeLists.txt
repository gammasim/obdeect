cmake_minimum_required(VERSION 3.20)

# Single source of truth for the version: the package version. Wheels and sdists get it
# from setuptools-scm (scikit-build-core passes SKBUILD_PROJECT_VERSION). Other builds use
# the newest v* git tag, or 0.0.0 without one. project() gets the numeric part (X.Y.Z).
if(DEFINED SKBUILD_PROJECT_VERSION)
  set(GSTRICKD_PACKAGE_VERSION "${SKBUILD_PROJECT_VERSION}")
else()
  execute_process(COMMAND git describe --tags --match "v[0-9]*" --abbrev=0
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  OUTPUT_VARIABLE GSTRICKD_GIT_TAG OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_QUIET RESULT_VARIABLE GSTRICKD_GIT_RESULT)
  if(GSTRICKD_GIT_RESULT EQUAL 0 AND GSTRICKD_GIT_TAG MATCHES "^v?([0-9]+\\.[0-9]+\\.[0-9]+)")
    set(GSTRICKD_PACKAGE_VERSION "${CMAKE_MATCH_1}")
  else()
    set(GSTRICKD_PACKAGE_VERSION "0.0.0")
  endif()
endif()
string(REGEX MATCH "^[0-9]+\\.[0-9]+\\.[0-9]+" GSTRICKD_NUMERIC_VERSION "${GSTRICKD_PACKAGE_VERSION}")
if(NOT GSTRICKD_NUMERIC_VERSION)
  message(FATAL_ERROR "gstrickd: cannot parse the package version '${GSTRICKD_PACKAGE_VERSION}'")
endif()
project(gstrickd VERSION ${GSTRICKD_NUMERIC_VERSION} LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(GSTRICKD_BUILD_TESTING "Build the C++ tests" ON)
option(GSTRICKD_BUILD_TOOLS "Build the command-line tools (ps-cat, ps-tee, corsika2ps, ps2simtel)" ON)
option(GSTRICKD_BUILD_CAPI "Build the C API shared library libgstrickd_c" ON)
option(GSTRICKD_BUILD_BENCH "Build the photon-stream benchmark" OFF)
if(DEFINED SKBUILD_PLATLIB_DIR)
  set(GSTRICKD_BUILD_PYTHON_API_DEFAULT ON)
else()
  set(GSTRICKD_BUILD_PYTHON_API_DEFAULT OFF)
endif()
option(GSTRICKD_BUILD_PYTHON_API "Build the Python extension" ${GSTRICKD_BUILD_PYTHON_API_DEFAULT})

# Install layout. In a wheel the headers and the CMake package live inside the Python
# package (python -m gstrickd --include-dir / --cmake-dir); otherwise they follow GNUInstallDirs.
if(DEFINED SKBUILD_PLATLIB_DIR)
  set(GSTRICKD_INCLUDE_INSTALL_DIR gstrickd/include)
  set(GSTRICKD_LIB_INSTALL_DIR gstrickd/lib)
  set(GSTRICKD_CMAKE_INSTALL_DIR gstrickd/cmake)
else()
  set(GSTRICKD_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR})
  set(GSTRICKD_LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})
  set(GSTRICKD_CMAKE_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/gstrickd)
endif()

# Generated headers (version.hpp) live in the build tree.
set(GSTRICKD_GENERATED_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.hpp.in
               ${GSTRICKD_GENERATED_DIR}/gstrickd/version.hpp @ONLY)

# Header-only, standard-library-only core. Exported as gstrickd::gstrickd_core; the in-tree
# name gstrickd_core stays valid.
add_library(gstrickd_core INTERFACE)
add_library(gstrickd::gstrickd_core ALIAS gstrickd_core)
target_include_directories(gstrickd_core INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cpp/include>
  $<BUILD_INTERFACE:${GSTRICKD_GENERATED_DIR}>
  $<INSTALL_INTERFACE:${GSTRICKD_INCLUDE_INSTALL_DIR}>)
target_compile_features(gstrickd_core INTERFACE cxx_std_20)
# Identical floating-point results for any thread count and chunk size on one
# platform: no implicit fused multiply-add contraction.
target_compile_options(gstrickd_core INTERFACE -ffp-contract=off)

function(gstrickd_enable_warnings target)
  target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
endfunction()

# Headers and the package config for consumers (find_package(gstrickd CONFIG)).
install(DIRECTORY cpp/include/gstrickd DESTINATION ${GSTRICKD_INCLUDE_INSTALL_DIR})
install(FILES ${GSTRICKD_GENERATED_DIR}/gstrickd/version.hpp
        DESTINATION ${GSTRICKD_INCLUDE_INSTALL_DIR}/gstrickd)
install(TARGETS gstrickd_core EXPORT gstrickd-targets)

# C API: stable ABI for non-C++ programs (cpp/include/gstrickd/gstrickd.h), exported as
# gstrickd::gstrickd_c.
if(GSTRICKD_BUILD_CAPI)
  add_library(gstrickd_c SHARED cpp/capi/gstrickd_c.cpp)
  add_library(gstrickd::gstrickd_c ALIAS gstrickd_c)
  target_link_libraries(gstrickd_c PRIVATE gstrickd_core)
  target_include_directories(gstrickd_c PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cpp/include>
    $<INSTALL_INTERFACE:${GSTRICKD_INCLUDE_INSTALL_DIR}>)
  target_compile_definitions(gstrickd_c PRIVATE GSTRICKD_C_BUILDING
                             "GSTRICKD_PACKAGE_VERSION=\"${GSTRICKD_PACKAGE_VERSION}\"")
  set_target_properties(gstrickd_c PROPERTIES CXX_VISIBILITY_PRESET hidden
                                              VISIBILITY_INLINES_HIDDEN ON)
  if(NOT DEFINED SKBUILD_PLATLIB_DIR) # wheels cannot hold the symlinks
    set_target_properties(gstrickd_c PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION 1)
  endif()
  gstrickd_enable_warnings(gstrickd_c)
  install(TARGETS gstrickd_c EXPORT gstrickd-targets
          LIBRARY DESTINATION ${GSTRICKD_LIB_INSTALL_DIR}
          RUNTIME DESTINATION ${GSTRICKD_LIB_INSTALL_DIR}
          ARCHIVE DESTINATION ${GSTRICKD_LIB_INSTALL_DIR})
endif()
install(EXPORT gstrickd-targets NAMESPACE gstrickd:: DESTINATION ${GSTRICKD_CMAKE_INSTALL_DIR})
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/gstrickd-config.cmake.in
                              ${CMAKE_CURRENT_BINARY_DIR}/gstrickd-config.cmake
                              INSTALL_DESTINATION ${GSTRICKD_CMAKE_INSTALL_DIR})
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/gstrickd-config-version.cmake
                                 VERSION ${PROJECT_VERSION} COMPATIBILITY SameMinorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gstrickd-config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/gstrickd-config-version.cmake
        DESTINATION ${GSTRICKD_CMAKE_INSTALL_DIR})

if(GSTRICKD_BUILD_PYTHON_API)
  find_package(Python 3.14 COMPONENTS Interpreter Development.Module REQUIRED)
  execute_process(COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
                  OUTPUT_VARIABLE nanobind_DIR OUTPUT_STRIP_TRAILING_WHITESPACE
                  COMMAND_ERROR_IS_FATAL ANY)
  find_package(nanobind CONFIG REQUIRED)
  nanobind_add_module(_core NB_STATIC cpp/bindings/core.cpp)
  target_link_libraries(_core PRIVATE gstrickd_core)
  gstrickd_enable_warnings(_core)
  install(TARGETS _core LIBRARY DESTINATION gstrickd)
endif()

if(GSTRICKD_BUILD_TOOLS)
  set(GSTRICKD_TOOLS ps-cat ps-tee corsika2ps ps2simtel)
  foreach(tool ${GSTRICKD_TOOLS})
    string(REPLACE "-" "_" source ${tool})
    add_executable(${tool} cpp/tools/${source}.cpp)
    target_link_libraries(${tool} PRIVATE gstrickd_core)
    target_compile_definitions(${tool} PRIVATE
                               "GSTRICKD_PACKAGE_VERSION=\"${GSTRICKD_PACKAGE_VERSION}\"")
    gstrickd_enable_warnings(${tool})
  endforeach()
  if(DEFINED SKBUILD_SCRIPTS_DIR)
    install(TARGETS ${GSTRICKD_TOOLS} RUNTIME DESTINATION ${SKBUILD_SCRIPTS_DIR})
  else()
    install(TARGETS ${GSTRICKD_TOOLS} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  endif()
endif()

if(GSTRICKD_BUILD_BENCH)
  foreach(bench bench_stream gen_corsika)
    add_executable(${bench} cpp/bench/${bench}.cpp)
    target_link_libraries(${bench} PRIVATE gstrickd_core)
    gstrickd_enable_warnings(${bench})
  endforeach()
endif()

if(GSTRICKD_BUILD_TESTING)
  enable_testing()
  if(GSTRICKD_BUILD_BENCH)
    add_test(NAME bench_stream_partial_chunk COMMAND bench_stream 65537
             ${CMAKE_CURRENT_BINARY_DIR}/bench_stream_test.ps)
  endif()
  if(GSTRICKD_BUILD_CAPI)
    enable_language(C)
    set(CMAKE_C_STANDARD 11)
    set(CMAKE_C_STANDARD_REQUIRED ON)
    set(CMAKE_C_EXTENSIONS OFF)
    add_executable(test_capi cpp/tests/test_capi.c)
    target_link_libraries(test_capi PRIVATE gstrickd_c)
    target_compile_options(test_capi PRIVATE -Wall -Wextra -Wpedantic -Werror)
    set(GSTRICKD_CONFORMANCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cpp/tests/data/conformance)
    add_test(NAME capi COMMAND test_capi ${CMAKE_CURRENT_BINARY_DIR}/capi_roundtrip.ps
             ${GSTRICKD_CONFORMANCE_DIR})
    set_tests_properties(capi PROPERTIES FIXTURES_SETUP capi_stream)
    add_executable(test_capi_cpp cpp/tests/test_capi_cpp.cpp)
    target_link_libraries(test_capi_cpp PRIVATE gstrickd_c gstrickd_core)
    gstrickd_enable_warnings(test_capi_cpp)
    add_test(NAME capi_cpp COMMAND test_capi_cpp ${CMAKE_CURRENT_BINARY_DIR}/capi_roundtrip.ps
             ${GSTRICKD_CONFORMANCE_DIR})
    set_tests_properties(capi_cpp PROPERTIES FIXTURES_REQUIRED capi_stream)
  endif()
  foreach(test_name version format stream eventio corsika simtel stages)
    add_executable(test_${test_name} cpp/tests/test_${test_name}.cpp)
    target_link_libraries(test_${test_name} PRIVATE gstrickd_core)
    gstrickd_enable_warnings(test_${test_name})
    add_test(NAME ${test_name} COMMAND test_${test_name})
  endforeach()

  # --version of each tool: "<tool> <package version> (photon-stream format <major>.<minor>)".
  file(READ ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include/gstrickd/format.hpp GSTRICKD_FORMAT_HPP)
  string(REGEX MATCH "format_major = ([0-9]+)" _ "${GSTRICKD_FORMAT_HPP}")
  set(GSTRICKD_FORMAT_MAJOR ${CMAKE_MATCH_1})
  string(REGEX MATCH "format_minor = ([0-9]+)" _ "${GSTRICKD_FORMAT_HPP}")
  set(GSTRICKD_FORMAT_MINOR ${CMAKE_MATCH_1})
  if(GSTRICKD_BUILD_TOOLS)
    foreach(tool ${GSTRICKD_TOOLS})
      add_test(NAME ${tool}_version COMMAND ${tool} --version)
      set_tests_properties(${tool}_version PROPERTIES PASS_REGULAR_EXPRESSION
        "^${tool} ${GSTRICKD_PACKAGE_VERSION} \\(photon-stream format ${GSTRICKD_FORMAT_MAJOR}\\.${GSTRICKD_FORMAT_MINOR}\\)")
    endforeach()
  endif()
endif()
