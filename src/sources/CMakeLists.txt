cmake_minimum_required(VERSION 3.20)
project(dlamb VERSION 0.1.0 LANGUAGES CXX)

option(DLAMB_BUILD_TESTING "Build the C++ tests" ON)
if(DEFINED SKBUILD_PLATLIB_DIR)
  set(DLAMB_BUILD_PYTHON_API_DEFAULT ON)
else()
  set(DLAMB_BUILD_PYTHON_API_DEFAULT OFF)
endif()
option(DLAMB_BUILD_PYTHON_API "Build the Python extension" ${DLAMB_BUILD_PYTHON_API_DEFAULT})

include(FetchContent)

# gstrickd (photon-stream format, Philox draws, stage ids). Lookup order:
#   1. an installed package: gstrickd_DIR or CMAKE_PREFIX_PATH;
#   2. the package of the Python interpreter: python -m gstrickd --cmake-dir;
#   3. FetchContent of a pinned commit (not with DLAMB_REQUIRE_INSTALLED_GSTRICKD=ON).
option(DLAMB_REQUIRE_INSTALLED_GSTRICKD
       "Require an installed gstrickd package (no FetchContent fallback)" OFF)
set(DLAMB_GSTRICKD_GIT_TAG 5af00790005283f81a46a4a49cfef97cd41556ff)

find_package(gstrickd 0.1 CONFIG QUIET)
if(NOT gstrickd_FOUND AND NOT gstrickd_DIR)
  find_package(Python 3.14 COMPONENTS Interpreter QUIET)
  if(Python_EXECUTABLE)
    execute_process(COMMAND "${Python_EXECUTABLE}" -m gstrickd --cmake-dir
                    OUTPUT_VARIABLE dlamb_gstrickd_python_dir
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    RESULT_VARIABLE dlamb_gstrickd_python_result
                    ERROR_QUIET)
    if(dlamb_gstrickd_python_result EQUAL 0 AND
       EXISTS "${dlamb_gstrickd_python_dir}/gstrickd-config.cmake")
      set(gstrickd_DIR "${dlamb_gstrickd_python_dir}")
      find_package(gstrickd 0.1 CONFIG QUIET)
    endif()
  endif()
endif()

if(gstrickd_FOUND)
  message(STATUS "gstrickd: installed package ${gstrickd_DIR} (version ${gstrickd_VERSION})")
else()
  if(DLAMB_REQUIRE_INSTALLED_GSTRICKD)
    message(FATAL_ERROR "gstrickd: no installed package found. Set gstrickd_DIR to the "
                        "directory of gstrickd-config.cmake, or install gstrickd.")
  endif()
  message(STATUS "gstrickd: no installed package found, fetching commit "
                 "${DLAMB_GSTRICKD_GIT_TAG} from GitHub (headers only)")
  # Only the source is fetched and not added as a subdirectory: the install rules, tests
  # and tools of the fetched project then stay out of dlamb's build and install.
  FetchContent_Populate(gstrickd
                        GIT_REPOSITORY https://github.com/GernotMaier/gstrickd.git
                        GIT_TAG ${DLAMB_GSTRICKD_GIT_TAG})
  add_library(gstrickd_fetched INTERFACE)
  target_include_directories(gstrickd_fetched INTERFACE
                             ${gstrickd_SOURCE_DIR}/cpp/include)
  target_compile_features(gstrickd_fetched INTERFACE cxx_std_20)
  add_library(gstrickd::gstrickd_core ALIAS gstrickd_fetched)
endif()

# Header-only, standard-library-only core. Photon chunks use the gstrickd headers.
add_library(dlamb_core INTERFACE)
target_include_directories(dlamb_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include)
target_compile_features(dlamb_core INTERFACE cxx_std_20)
target_link_libraries(dlamb_core INTERFACE gstrickd::gstrickd_core)
# Identical floating-point results for any thread count and chunk size on one
# platform: no implicit fused multiply-add contraction.
target_compile_options(dlamb_core INTERFACE -ffp-contract=off)

function(dlamb_enable_warnings target)
  target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
  # GCC 13 reports false maybe-uninitialized warnings inside std::variant at -O3.
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(${target} PRIVATE -Wno-error=maybe-uninitialized)
  endif()
endfunction()

if(DLAMB_BUILD_PYTHON_API)
  find_package(Python 3.14 COMPONENTS Interpreter Development.Module REQUIRED)
  execute_process(COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
                  OUTPUT_VARIABLE nanobind_DIR OUTPUT_STRIP_TRAILING_WHITESPACE
                  COMMAND_ERROR_IS_FATAL ANY)
  find_package(nanobind CONFIG REQUIRED)
  nanobind_add_module(_core NB_STATIC cpp/bindings/core.cpp)
  target_link_libraries(_core PRIVATE dlamb_core)
  dlamb_enable_warnings(_core)
  install(TARGETS _core LIBRARY DESTINATION dlamb)
endif()

if(DLAMB_BUILD_TESTING)
  enable_testing()
  foreach(test_name version density source gstrickd_contract)
    add_executable(test_${test_name} cpp/tests/test_${test_name}.cpp)
    target_link_libraries(test_${test_name} PRIVATE dlamb_core)
    dlamb_enable_warnings(test_${test_name})
    add_test(NAME ${test_name} COMMAND test_${test_name})
  endforeach()
endif()

# Light-source programs: tool and test sources share the program name.
option(DLAMB_BUILD_TOOLS "Build the light-source command-line tools" ON)
set(DLAMB_LIGHT_SOURCE_PROGRAMS
    flasher led point-source beam-scan laser-beam night-sky muon-track)
if(DLAMB_BUILD_TESTING)
  foreach(test_name ${DLAMB_LIGHT_SOURCE_PROGRAMS} sky_args)
    add_executable(test_${test_name} cpp/tests/test_${test_name}.cpp)
    target_link_libraries(test_${test_name} PRIVATE dlamb_core)
    dlamb_enable_warnings(test_${test_name})
    add_test(NAME ${test_name} COMMAND test_${test_name})
  endforeach()
endif()
if(DLAMB_BUILD_TOOLS)
  foreach(tool ${DLAMB_LIGHT_SOURCE_PROGRAMS})
    add_executable(dlamb-${tool} cpp/tools/${tool}.cpp)
    target_link_libraries(dlamb-${tool} PRIVATE dlamb_core)
    dlamb_enable_warnings(dlamb-${tool})
  endforeach()
  if(DEFINED SKBUILD_SCRIPTS_DIR)
    foreach(tool ${DLAMB_LIGHT_SOURCE_PROGRAMS})
      install(TARGETS dlamb-${tool} RUNTIME DESTINATION ${SKBUILD_SCRIPTS_DIR})
    endforeach()
  endif()
endif()
