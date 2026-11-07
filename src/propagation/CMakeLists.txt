cmake_minimum_required(VERSION 3.20)
project(dlufd VERSION 0.1.0 LANGUAGES CXX)

option(DLUFD_BUILD_TESTING "Build the C++ tests" ON)
if(DEFINED SKBUILD_PLATLIB_DIR)
  set(DLUFD_BUILD_PYTHON_API_DEFAULT ON)
else()
  set(DLUFD_BUILD_PYTHON_API_DEFAULT OFF)
endif()
option(DLUFD_BUILD_PYTHON_API "Build the Python extension" ${DLUFD_BUILD_PYTHON_API_DEFAULT})

include(FetchContent)
FetchContent_Declare(gstrickd GIT_REPOSITORY https://github.com/GernotMaier/gstrickd.git
                     GIT_TAG fcffa70530510ad013edf644f700ed9989912682)
set(GSTRICKD_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(GSTRICKD_BUILD_PYTHON_API OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(gstrickd)

# Header-only, standard-library-only core.
add_library(dlufd_core INTERFACE)
target_include_directories(dlufd_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include)
target_compile_features(dlufd_core INTERFACE cxx_std_20)
target_link_libraries(dlufd_core INTERFACE gstrickd_core)
# Identical floating-point results for any thread count and chunk size on one
# platform: no implicit fused multiply-add contraction.
if(NOT MSVC)
  target_compile_options(dlufd_core INTERFACE -ffp-contract=off)
endif()

function(dlufd_enable_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /WX)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()
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

if(DLUFD_BUILD_TESTING)
  enable_testing()
  foreach(test_name version transmission thinning pipeline)
    add_executable(test_${test_name} cpp/tests/test_${test_name}.cpp)
    target_link_libraries(test_${test_name} PRIVATE dlufd_core)
    # Reference data: cpp/tests/data (copied from sim_telarray, see test_transmission.cpp).
    target_compile_definitions(test_${test_name}
                               PRIVATE DLUFD_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/cpp/tests/data")
    dlufd_enable_warnings(test_${test_name})
    add_test(NAME ${test_name} COMMAND test_${test_name})
  endforeach()
endif()
