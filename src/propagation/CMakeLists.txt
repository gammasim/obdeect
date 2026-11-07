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

# gstrickd: an installed package first (find_package, then the gstrickd of the Python
# environment that runs CMake), FetchContent at a pinned tag only as a fallback.
option(DLUFD_REQUIRE_INSTALLED_GSTRICKD
       "Require an installed gstrickd package and never fetch it (offline and sdist builds)" OFF)
find_package(gstrickd 0.1 CONFIG QUIET)
if(NOT gstrickd_FOUND AND NOT gstrickd_DIR)
  find_package(Python 3.14 COMPONENTS Interpreter)
  if(Python_Interpreter_FOUND)
    execute_process(COMMAND "${Python_EXECUTABLE}" -m gstrickd --cmake-dir
                    OUTPUT_VARIABLE dlufd_gstrickd_cmake_dir OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET RESULT_VARIABLE dlufd_gstrickd_probe)
    if(dlufd_gstrickd_probe EQUAL 0 AND EXISTS "${dlufd_gstrickd_cmake_dir}/gstrickd-config.cmake")
      set(gstrickd_DIR "${dlufd_gstrickd_cmake_dir}")
      find_package(gstrickd 0.1 CONFIG QUIET)
    endif()
  endif()
endif()

if(gstrickd_FOUND)
  message(STATUS "dlufd: using installed gstrickd ${gstrickd_VERSION} (${gstrickd_DIR})")
else()
  if(DLUFD_REQUIRE_INSTALLED_GSTRICKD)
    message(FATAL_ERROR "dlufd: no installed gstrickd package. Set -Dgstrickd_DIR=<directory of "
                        "gstrickd-config.cmake> (python -m gstrickd --cmake-dir) or pip install gstrickd.")
  endif()
  message(STATUS "dlufd: no installed gstrickd package, fetching gstrickd 5af0079 (FetchContent)")
  include(FetchContent)
  # Headers only, fetched directly and not added as a subdirectory: gstrickd's install
  # rules, tests and tools then stay out of dlufd's build, wheel and installs.
  FetchContent_Populate(gstrickd GIT_REPOSITORY https://github.com/GernotMaier/gstrickd.git
                        GIT_TAG 5af00790005283f81a46a4a49cfef97cd41556ff)
  add_library(gstrickd_fetched INTERFACE)
  target_include_directories(gstrickd_fetched INTERFACE "${gstrickd_SOURCE_DIR}/cpp/include")
  target_compile_features(gstrickd_fetched INTERFACE cxx_std_20)
  add_library(gstrickd::gstrickd_core ALIAS gstrickd_fetched)
endif()

# Header-only, standard-library-only core.
add_library(dlufd_core INTERFACE)
target_include_directories(dlufd_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include)
target_compile_features(dlufd_core INTERFACE cxx_std_20)
target_link_libraries(dlufd_core INTERFACE gstrickd::gstrickd_core)
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
                    propagator scattering gstrickd_contract)
    add_executable(test_${test_name} cpp/tests/test_${test_name}.cpp)
    target_link_libraries(test_${test_name} PRIVATE dlufd_core)
    # Reference data: cpp/tests/data (copied from sim_telarray, see test_transmission.cpp).
    target_compile_definitions(test_${test_name}
                               PRIVATE DLUFD_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/cpp/tests/data")
    dlufd_enable_warnings(test_${test_name})
    add_test(NAME ${test_name} COMMAND test_${test_name})
  endforeach()

  # Agreement with CORSIKA's atmo.c: DLUFD_CORSIKA_REFERENCE=<corsika7/bernlohr directory>.
  option(DLUFD_CORSIKA_REFERENCE "Path to CORSIKA's corsika7/bernlohr directory (atmo.c)" OFF)
  if(DLUFD_CORSIKA_REFERENCE)
    enable_language(C)
    add_library(corsika_atmo STATIC
      ${DLUFD_CORSIKA_REFERENCE}/atmo.c
      ${DLUFD_CORSIKA_REFERENCE}/fileopen.c
      ${DLUFD_CORSIKA_REFERENCE}/straux.c)
    target_include_directories(corsika_atmo PUBLIC ${DLUFD_CORSIKA_REFERENCE})
    target_compile_options(corsika_atmo PRIVATE -w -ffp-contract=off)
    add_executable(test_corsika_reference cpp/tests/test_corsika_reference.cpp)
    target_link_libraries(test_corsika_reference PRIVATE dlufd_core corsika_atmo m)
    target_compile_definitions(test_corsika_reference PRIVATE
      DLUFD_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/cpp/tests/data")
    add_test(NAME corsika_reference COMMAND test_corsika_reference)
  endif()
endif()
