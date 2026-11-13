cmake_minimum_required(VERSION 3.20)
project(obdeect_sources LANGUAGES CXX)

option(OBDEECT_SOURCES_BUILD_TESTING "Build the C++ tests" ON)
if(DEFINED SKBUILD_PLATLIB_DIR)
  set(OBDEECT_SOURCES_BUILD_PYTHON_API_DEFAULT ON)
else()
  set(OBDEECT_SOURCES_BUILD_PYTHON_API_DEFAULT OFF)
endif()
option(OBDEECT_SOURCES_BUILD_PYTHON_API "Build the Python extension" ${OBDEECT_SOURCES_BUILD_PYTHON_API_DEFAULT})

# The stream component: the target of the combined build, else the stream sources next to
# this component (builds from the repository), else an installed package (sdist builds).
if(NOT TARGET obdeect::stream)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../stream/CMakeLists.txt)
    set(OBDEECT_STREAM_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(OBDEECT_STREAM_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
    set(OBDEECT_STREAM_BUILD_CAPI OFF CACHE BOOL "" FORCE)
    set(OBDEECT_STREAM_BUILD_PYTHON_API OFF CACHE BOOL "" FORCE)
    set(OBDEECT_STREAM_INSTALL OFF CACHE BOOL "" FORCE)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../stream ${CMAKE_CURRENT_BINARY_DIR}/stream)
    message(STATUS "obdeect sources: using the stream sources in ../stream")
  else()
    find_package(Python 3.14 COMPONENTS Interpreter)
    if(Python_Interpreter_FOUND AND NOT obdeect_stream_DIR)
      execute_process(COMMAND "${Python_EXECUTABLE}" -m obdeect.stream --cmake-dir
                      OUTPUT_VARIABLE obdeect_stream_DIR OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    endif()
    find_package(obdeect_stream CONFIG REQUIRED)
    message(STATUS "obdeect sources: using the installed stream package ${obdeect_stream_VERSION} (${obdeect_stream_DIR})")
  endif()
endif()

# Header-only, standard-library-only core. Photon chunks use the obdeect stream headers.
add_library(obdeect_sources_core INTERFACE)
target_include_directories(obdeect_sources_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include)
target_compile_features(obdeect_sources_core INTERFACE cxx_std_20)
target_link_libraries(obdeect_sources_core INTERFACE obdeect::stream)
# Identical floating-point results for any thread count and chunk size on one
# platform: no implicit fused multiply-add contraction.
target_compile_options(obdeect_sources_core INTERFACE -ffp-contract=off)

function(obdeect_sources_enable_warnings target)
  target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
  # GCC 13 reports false maybe-uninitialized warnings inside std::variant at -O3.
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(${target} PRIVATE -Wno-error=maybe-uninitialized)
  endif()
endfunction()

if(OBDEECT_SOURCES_BUILD_PYTHON_API)
  find_package(Python 3.14 COMPONENTS Interpreter Development.Module REQUIRED)
  execute_process(COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
                  OUTPUT_VARIABLE nanobind_DIR OUTPUT_STRIP_TRAILING_WHITESPACE
                  COMMAND_ERROR_IS_FATAL ANY)
  find_package(nanobind CONFIG REQUIRED)
  nanobind_add_module(_core NB_STATIC cpp/bindings/core.cpp)
  target_link_libraries(_core PRIVATE obdeect_sources_core)
  obdeect_sources_enable_warnings(_core)
  install(TARGETS _core LIBRARY DESTINATION obdeect/sources)
endif()

if(OBDEECT_SOURCES_BUILD_TESTING)
  enable_testing()
  foreach(test_name version density source stream_contract)
    add_executable(obdeect_sources_test_${test_name} cpp/tests/test_${test_name}.cpp)
    target_link_libraries(obdeect_sources_test_${test_name} PRIVATE obdeect_sources_core)
    # check.hpp lives with the stream tests; tests are built from the repository only, never from an sdist.
    target_include_directories(obdeect_sources_test_${test_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../stream/cpp/tests)
    obdeect_sources_enable_warnings(obdeect_sources_test_${test_name})
    add_test(NAME ${test_name} COMMAND obdeect_sources_test_${test_name})
  endforeach()
endif()

# Light-source programs: tool and test sources share the program name.
option(OBDEECT_SOURCES_BUILD_TOOLS "Build the light-source command-line tools" ON)
set(OBDEECT_SOURCES_LIGHT_SOURCE_PROGRAMS
    flasher led point-source beam-scan laser-beam night-sky muon-track)
if(OBDEECT_SOURCES_BUILD_TESTING)
  foreach(test_name ${OBDEECT_SOURCES_LIGHT_SOURCE_PROGRAMS} sky_args)
    add_executable(obdeect_sources_test_${test_name} cpp/tests/test_${test_name}.cpp)
    target_link_libraries(obdeect_sources_test_${test_name} PRIVATE obdeect_sources_core)
    # check.hpp lives with the stream tests; tests are built from the repository only, never from an sdist.
    target_include_directories(obdeect_sources_test_${test_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../stream/cpp/tests)
    obdeect_sources_enable_warnings(obdeect_sources_test_${test_name})
    add_test(NAME ${test_name} COMMAND obdeect_sources_test_${test_name})
  endforeach()
endif()
if(OBDEECT_SOURCES_BUILD_TOOLS)
  foreach(tool ${OBDEECT_SOURCES_LIGHT_SOURCE_PROGRAMS})
    add_executable(obdeect-sources-${tool} cpp/tools/${tool}.cpp)
    target_link_libraries(obdeect-sources-${tool} PRIVATE obdeect_sources_core)
    obdeect_sources_enable_warnings(obdeect-sources-${tool})
  endforeach()
  # Wheels: inside the package, run through `obdeect sources <program>`.
  if(DEFINED SKBUILD_PLATLIB_DIR)
    foreach(tool ${OBDEECT_SOURCES_LIGHT_SOURCE_PROGRAMS})
      install(TARGETS obdeect-sources-${tool} RUNTIME DESTINATION obdeect/sources/bin)
    endforeach()
  endif()
endif()
