cmake_minimum_required(VERSION 3.20)
project(gstrickd VERSION 0.1.0 LANGUAGES CXX)

option(GSTRICKD_BUILD_TESTING "Build the C++ tests" ON)
if(DEFINED SKBUILD_PLATLIB_DIR)
  set(GSTRICKD_BUILD_PYTHON_API_DEFAULT ON)
else()
  set(GSTRICKD_BUILD_PYTHON_API_DEFAULT OFF)
endif()
option(GSTRICKD_BUILD_PYTHON_API "Build the Python extension" ${GSTRICKD_BUILD_PYTHON_API_DEFAULT})

# Header-only, standard-library-only core.
add_library(gstrickd_core INTERFACE)
target_include_directories(gstrickd_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include)
target_compile_features(gstrickd_core INTERFACE cxx_std_20)
# Identical floating-point results for any thread count and chunk size on one
# platform: no implicit fused multiply-add contraction.
if(NOT MSVC)
  target_compile_options(gstrickd_core INTERFACE -ffp-contract=off)
endif()

function(gstrickd_enable_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /WX)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()
endfunction()

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

if(GSTRICKD_BUILD_TESTING)
  enable_testing()
  foreach(test_name version)
    add_executable(test_${test_name} cpp/tests/test_${test_name}.cpp)
    target_link_libraries(test_${test_name} PRIVATE gstrickd_core)
    gstrickd_enable_warnings(test_${test_name})
    add_test(NAME ${test_name} COMMAND test_${test_name})
  endforeach()
endif()
