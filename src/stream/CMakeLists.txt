cmake_minimum_required(VERSION 3.20)

# Single source of truth for the version: the package version. Wheels and sdists get it
# from setuptools-scm (scikit-build-core passes SKBUILD_PROJECT_VERSION). Other builds use
# the newest v* git tag, or 0.0.0 without one. project() gets the numeric part (X.Y.Z).
if(DEFINED SKBUILD_PROJECT_VERSION)
  set(OBDEECT_STREAM_PACKAGE_VERSION "${SKBUILD_PROJECT_VERSION}")
else()
  execute_process(COMMAND git describe --tags --match "v[0-9]*" --abbrev=0
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  OUTPUT_VARIABLE OBDEECT_STREAM_GIT_TAG OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_QUIET RESULT_VARIABLE OBDEECT_STREAM_GIT_RESULT)
  if(OBDEECT_STREAM_GIT_RESULT EQUAL 0 AND OBDEECT_STREAM_GIT_TAG MATCHES "^v?([0-9]+\\.[0-9]+\\.[0-9]+)")
    set(OBDEECT_STREAM_PACKAGE_VERSION "${CMAKE_MATCH_1}")
  else()
    set(OBDEECT_STREAM_PACKAGE_VERSION "0.0.0")
  endif()
endif()
string(REGEX MATCH "^[0-9]+\\.[0-9]+\\.[0-9]+" OBDEECT_STREAM_NUMERIC_VERSION "${OBDEECT_STREAM_PACKAGE_VERSION}")
if(NOT OBDEECT_STREAM_NUMERIC_VERSION)
  message(FATAL_ERROR "obdeect stream: cannot parse the package version '${OBDEECT_STREAM_PACKAGE_VERSION}'")
endif()
project(obdeect_stream VERSION ${OBDEECT_STREAM_NUMERIC_VERSION} LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(OBDEECT_STREAM_BUILD_TESTING "Build the C++ tests" ON)
option(OBDEECT_STREAM_BUILD_TOOLS "Build the command-line tools of obdeect stream (cat, tee, from-corsika, to-simtel, to-grisu)" ON)
option(OBDEECT_STREAM_BUILD_CAPI "Build the C API shared library libobdeect_stream_c" ON)
option(OBDEECT_STREAM_INSTALL "Install the headers and the CMake package (OFF when another component embeds the stream sources)" ON)
option(OBDEECT_STREAM_BUILD_BENCH "Build the photon-stream benchmark" OFF)
if(DEFINED SKBUILD_PLATLIB_DIR)
  set(OBDEECT_STREAM_BUILD_PYTHON_API_DEFAULT ON)
else()
  set(OBDEECT_STREAM_BUILD_PYTHON_API_DEFAULT OFF)
endif()
option(OBDEECT_STREAM_BUILD_PYTHON_API "Build the Python extension" ${OBDEECT_STREAM_BUILD_PYTHON_API_DEFAULT})

# Install layout. In a wheel the headers and the CMake package live inside the Python
# package (python -m obdeect.stream --include-dir / --cmake-dir); otherwise they follow GNUInstallDirs.
if(DEFINED SKBUILD_PLATLIB_DIR)
  set(OBDEECT_STREAM_INCLUDE_INSTALL_DIR obdeect/stream/include)
  set(OBDEECT_STREAM_LIB_INSTALL_DIR obdeect/stream/lib)
  set(OBDEECT_STREAM_CMAKE_INSTALL_DIR obdeect/stream/cmake)
else()
  set(OBDEECT_STREAM_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR})
  set(OBDEECT_STREAM_LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})
  set(OBDEECT_STREAM_CMAKE_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/obdeect_stream)
endif()

# Generated headers (version.hpp) live in the build tree.
set(OBDEECT_STREAM_GENERATED_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.hpp.in
               ${OBDEECT_STREAM_GENERATED_DIR}/obdeect/stream/version.hpp @ONLY)

# Header-only, standard-library-only core. Exported as obdeect::stream; the in-tree
# name obdeect_stream_core stays valid.
add_library(obdeect_stream_core INTERFACE)
add_library(obdeect::stream ALIAS obdeect_stream_core)
set_target_properties(obdeect_stream_core PROPERTIES EXPORT_NAME stream)
target_include_directories(obdeect_stream_core INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cpp/include>
  $<BUILD_INTERFACE:${OBDEECT_STREAM_GENERATED_DIR}>
  $<INSTALL_INTERFACE:${OBDEECT_STREAM_INCLUDE_INSTALL_DIR}>)
target_compile_features(obdeect_stream_core INTERFACE cxx_std_20)
# Identical floating-point results for any thread count and chunk size on one
# platform: no implicit fused multiply-add contraction.
target_compile_options(obdeect_stream_core INTERFACE -ffp-contract=off)

function(obdeect_stream_enable_warnings target)
  target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
endfunction()

# Headers and the package config for consumers (find_package(obdeect_stream CONFIG)).
if(OBDEECT_STREAM_INSTALL)
  install(DIRECTORY cpp/include/obdeect/stream DESTINATION ${OBDEECT_STREAM_INCLUDE_INSTALL_DIR}/obdeect)
  install(FILES ${OBDEECT_STREAM_GENERATED_DIR}/obdeect/stream/version.hpp
          DESTINATION ${OBDEECT_STREAM_INCLUDE_INSTALL_DIR}/obdeect/stream)
  install(TARGETS obdeect_stream_core EXPORT obdeect_stream-targets)
endif()

# C API: stable ABI for non-C++ programs (cpp/include/obdeect/stream/obdeect_stream.h), exported as
# obdeect::stream_c.
if(OBDEECT_STREAM_BUILD_CAPI)
  add_library(obdeect_stream_c SHARED cpp/capi/obdeect_stream_c.cpp)
  add_library(obdeect::stream_c ALIAS obdeect_stream_c)
  set_target_properties(obdeect_stream_c PROPERTIES EXPORT_NAME stream_c)
  target_link_libraries(obdeect_stream_c PRIVATE obdeect_stream_core)
  target_include_directories(obdeect_stream_c PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cpp/include>
    $<INSTALL_INTERFACE:${OBDEECT_STREAM_INCLUDE_INSTALL_DIR}>)
  target_compile_definitions(obdeect_stream_c PRIVATE OBDEECT_STREAM_C_BUILDING
                             "OBDEECT_STREAM_PACKAGE_VERSION=\"${OBDEECT_STREAM_PACKAGE_VERSION}\"")
  set_target_properties(obdeect_stream_c PROPERTIES CXX_VISIBILITY_PRESET hidden
                                              VISIBILITY_INLINES_HIDDEN ON)
  if(NOT DEFINED SKBUILD_PLATLIB_DIR) # wheels cannot hold the symlinks
    set_target_properties(obdeect_stream_c PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION 1)
  endif()
  obdeect_stream_enable_warnings(obdeect_stream_c)
  install(TARGETS obdeect_stream_c EXPORT obdeect_stream-targets
          LIBRARY DESTINATION ${OBDEECT_STREAM_LIB_INSTALL_DIR}
          RUNTIME DESTINATION ${OBDEECT_STREAM_LIB_INSTALL_DIR}
          ARCHIVE DESTINATION ${OBDEECT_STREAM_LIB_INSTALL_DIR})
endif()
if(OBDEECT_STREAM_INSTALL)
  install(EXPORT obdeect_stream-targets NAMESPACE obdeect:: DESTINATION ${OBDEECT_STREAM_CMAKE_INSTALL_DIR})
  configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/obdeect_stream-config.cmake.in
                                ${CMAKE_CURRENT_BINARY_DIR}/obdeect_stream-config.cmake
                                INSTALL_DESTINATION ${OBDEECT_STREAM_CMAKE_INSTALL_DIR})
  write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/obdeect_stream-config-version.cmake
                                   VERSION ${PROJECT_VERSION} COMPATIBILITY SameMinorVersion)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/obdeect_stream-config.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/obdeect_stream-config-version.cmake
          DESTINATION ${OBDEECT_STREAM_CMAKE_INSTALL_DIR})
endif()

if(OBDEECT_STREAM_BUILD_PYTHON_API)
  find_package(Python 3.14 COMPONENTS Interpreter Development.Module REQUIRED)
  execute_process(COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
                  OUTPUT_VARIABLE nanobind_DIR OUTPUT_STRIP_TRAILING_WHITESPACE
                  COMMAND_ERROR_IS_FATAL ANY)
  find_package(nanobind CONFIG REQUIRED)
  nanobind_add_module(_core NB_STATIC cpp/bindings/core.cpp)
  target_link_libraries(_core PRIVATE obdeect_stream_core)
  obdeect_stream_enable_warnings(_core)
  install(TARGETS _core LIBRARY DESTINATION obdeect/stream)
endif()

if(OBDEECT_STREAM_BUILD_TOOLS)
  # Commands of `obdeect stream <command>` and their sources in cpp/tools.
  set(OBDEECT_STREAM_COMMANDS cat tee from-corsika to-simtel to-grisu)
  set(OBDEECT_STREAM_COMMAND_SOURCES ps_cat ps_tee corsika2ps ps2simtel ps2grisu)
  set(OBDEECT_STREAM_TOOLS)
  foreach(command source IN ZIP_LISTS OBDEECT_STREAM_COMMANDS OBDEECT_STREAM_COMMAND_SOURCES)
    set(tool obdeect-stream-${command})
    list(APPEND OBDEECT_STREAM_TOOLS ${tool})
    add_executable(${tool} cpp/tools/${source}.cpp)
    target_link_libraries(${tool} PRIVATE obdeect_stream_core)
    target_compile_definitions(${tool} PRIVATE
                               "OBDEECT_STREAM_PACKAGE_VERSION=\"${OBDEECT_STREAM_PACKAGE_VERSION}\"")
    obdeect_stream_enable_warnings(${tool})
  endforeach()
  # Wheels: inside the package, run through `obdeect stream <command>`.
  if(DEFINED SKBUILD_PLATLIB_DIR)
    install(TARGETS ${OBDEECT_STREAM_TOOLS} RUNTIME DESTINATION obdeect/stream/bin)
  else()
    install(TARGETS ${OBDEECT_STREAM_TOOLS} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  endif()
endif()

if(OBDEECT_STREAM_BUILD_BENCH)
  foreach(bench bench_stream gen_corsika)
    add_executable(${bench} cpp/bench/${bench}.cpp)
    target_link_libraries(${bench} PRIVATE obdeect_stream_core)
    obdeect_stream_enable_warnings(${bench})
  endforeach()
endif()

if(OBDEECT_STREAM_BUILD_TESTING)
  enable_testing()
  if(OBDEECT_STREAM_BUILD_BENCH)
    add_test(NAME bench_stream_partial_chunk COMMAND bench_stream 65537
             ${CMAKE_CURRENT_BINARY_DIR}/bench_stream_test.pst)
  endif()
  if(OBDEECT_STREAM_BUILD_CAPI)
    enable_language(C)
    set(CMAKE_C_STANDARD 11)
    set(CMAKE_C_STANDARD_REQUIRED ON)
    set(CMAKE_C_EXTENSIONS OFF)
    add_executable(test_capi cpp/tests/test_capi.c)
    target_link_libraries(test_capi PRIVATE obdeect_stream_c)
    target_compile_options(test_capi PRIVATE -Wall -Wextra -Wpedantic -Werror)
    set(OBDEECT_STREAM_CONFORMANCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cpp/tests/data/conformance)
    add_test(NAME capi COMMAND test_capi ${CMAKE_CURRENT_BINARY_DIR}/capi_roundtrip.pst
             ${OBDEECT_STREAM_CONFORMANCE_DIR})
    set_tests_properties(capi PROPERTIES FIXTURES_SETUP capi_stream)
    add_executable(test_capi_cpp cpp/tests/test_capi_cpp.cpp)
    target_link_libraries(test_capi_cpp PRIVATE obdeect_stream_c obdeect_stream_core)
    obdeect_stream_enable_warnings(test_capi_cpp)
    add_test(NAME capi_cpp COMMAND test_capi_cpp ${CMAKE_CURRENT_BINARY_DIR}/capi_roundtrip.pst
             ${OBDEECT_STREAM_CONFORMANCE_DIR})
    set_tests_properties(capi_cpp PROPERTIES FIXTURES_REQUIRED capi_stream)
  endif()
  foreach(test_name version format parse cli stream eventio corsika simtel grisu stages wavelength)
    add_executable(test_${test_name} cpp/tests/test_${test_name}.cpp)
    target_link_libraries(test_${test_name} PRIVATE obdeect_stream_core)
    obdeect_stream_enable_warnings(test_${test_name})
    add_test(NAME ${test_name} COMMAND test_${test_name})
  endforeach()

  # --version of each tool: "<tool> <package version> (photon-stream format <major>.<minor>)".
  file(READ ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include/obdeect/stream/format.hpp OBDEECT_STREAM_FORMAT_HPP)
  string(REGEX MATCH "format_major = ([0-9]+)" _ "${OBDEECT_STREAM_FORMAT_HPP}")
  set(OBDEECT_STREAM_FORMAT_MAJOR ${CMAKE_MATCH_1})
  string(REGEX MATCH "format_minor = ([0-9]+)" _ "${OBDEECT_STREAM_FORMAT_HPP}")
  set(OBDEECT_STREAM_FORMAT_MINOR ${CMAKE_MATCH_1})
  if(OBDEECT_STREAM_BUILD_TOOLS)
    foreach(command ${OBDEECT_STREAM_COMMANDS})
      add_test(NAME ${command}_version COMMAND obdeect-stream-${command} --version)
      set_tests_properties(${command}_version PROPERTIES PASS_REGULAR_EXPRESSION
        "^obdeect stream ${command} ${OBDEECT_STREAM_PACKAGE_VERSION} \\(photon-stream format ${OBDEECT_STREAM_FORMAT_MAJOR}\\.${OBDEECT_STREAM_FORMAT_MINOR}\\)")
    endforeach()
  endif()
endif()
