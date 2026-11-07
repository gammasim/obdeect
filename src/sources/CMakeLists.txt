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
FetchContent_Declare(gstrickd GIT_REPOSITORY https://github.com/GernotMaier/gstrickd.git
                     GIT_TAG 9ac0c0ee7c8edce7f24d190cedda6d8c67397f19)
set(GSTRICKD_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(GSTRICKD_BUILD_PYTHON_API OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(gstrickd)

# Header-only, standard-library-only core. Photon chunks use the gstrickd headers.
add_library(dlamb_core INTERFACE)
target_include_directories(dlamb_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include)
target_compile_features(dlamb_core INTERFACE cxx_std_20)
target_link_libraries(dlamb_core INTERFACE gstrickd_core)
# Identical floating-point results for any thread count and chunk size on one
# platform: no implicit fused multiply-add contraction.
if(NOT MSVC)
  target_compile_options(dlamb_core INTERFACE -ffp-contract=off)
endif()

function(dlamb_enable_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /WX)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
    # GCC 13 reports false maybe-uninitialized warnings inside std::variant at -O3.
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      target_compile_options(${target} PRIVATE -Wno-error=maybe-uninitialized)
    endif()
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
  foreach(test_name version density source)
    add_executable(test_${test_name} cpp/tests/test_${test_name}.cpp)
    target_link_libraries(test_${test_name} PRIVATE dlamb_core)
    dlamb_enable_warnings(test_${test_name})
    add_test(NAME ${test_name} COMMAND test_${test_name})
  endforeach()
endif()

# LightEmission programs (tool and test sources named as the programs, or with
# underscores for fake-muon, octo, ls-test and ls-geo-test)
option(DLAMB_BUILD_TOOLS "Build the light-source command-line tools" ON)
if(DLAMB_BUILD_TESTING)
  foreach(test_name ff-1m ff-gct fpls pixled xyzls xyls ls-beam nsbls fake_muon octo ls_test
                    ls_geo_test)
    add_executable(test_${test_name} cpp/tests/test_${test_name}.cpp)
    target_link_libraries(test_${test_name} PRIVATE dlamb_core)
    dlamb_enable_warnings(test_${test_name})
    add_test(NAME ${test_name} COMMAND test_${test_name})
  endforeach()
endif()
if(DLAMB_BUILD_TOOLS)
  foreach(tool ff-1m ff-gct fpls pixled xyzls xyls ls-beam nsbls fake-muon octo ls-test
               ls-geo-test)
    string(REPLACE "-" "_" source ${tool})
    if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/cpp/tools/${source}.cpp)
      set(source ${tool})
    endif()
    add_executable(dlamb-${tool} cpp/tools/${source}.cpp)
    target_link_libraries(dlamb-${tool} PRIVATE dlamb_core)
    dlamb_enable_warnings(dlamb-${tool})
  endforeach()
  if(DEFINED SKBUILD_SCRIPTS_DIR)
    install(TARGETS dlamb-ff-1m dlamb-ff-gct dlamb-fpls dlamb-pixled dlamb-xyzls dlamb-xyls
                    dlamb-ls-beam dlamb-nsbls dlamb-fake-muon dlamb-octo dlamb-ls-test
                    dlamb-ls-geo-test
            RUNTIME DESTINATION ${SKBUILD_SCRIPTS_DIR})
  endif()
endif()
