cmake_minimum_required(VERSION 3.20)
project(obdeect_sources LANGUAGES CXX)

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

option(OBDEECT_SOURCES_BUILD_TESTING "Build the C++ tests" ON)
if(DEFINED SKBUILD_PLATLIB_DIR)
  set(OBDEECT_SOURCES_BUILD_PYTHON_API_DEFAULT ON)
else()
  set(OBDEECT_SOURCES_BUILD_PYTHON_API_DEFAULT OFF)
endif()
option(OBDEECT_SOURCES_BUILD_PYTHON_API "Build the Python extension" ${OBDEECT_SOURCES_BUILD_PYTHON_API_DEFAULT})

# The stream component (see obdeect_locate_stream in cmake/ObdeectCommon.cmake).
obdeect_locate_stream(sources)

# Header-only, standard-library-only core. Photon chunks use the obdeect stream headers.
add_library(obdeect_sources_core INTERFACE)
target_include_directories(obdeect_sources_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include)
target_link_libraries(obdeect_sources_core INTERFACE obdeect::stream)
obdeect_set_core_options(obdeect_sources_core)

function(obdeect_sources_enable_warnings target)
  obdeect_enable_warnings(${target})
  # GCC 13 reports false maybe-uninitialized warnings inside std::variant at -O3.
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(${target} PRIVATE -Wno-error=maybe-uninitialized)
  endif()
endfunction()

if(OBDEECT_SOURCES_BUILD_PYTHON_API)
  obdeect_add_python_module(sources obdeect_sources_core)
  obdeect_sources_enable_warnings(_core)
endif()

if(OBDEECT_SOURCES_BUILD_TESTING)
  enable_testing()
  foreach(test_name version density source emit run stream_contract)
    add_executable(obdeect_sources_test_${test_name} cpp/tests/test_${test_name}.cpp)
    target_link_libraries(obdeect_sources_test_${test_name} PRIVATE obdeect_sources_core)
    # check.hpp lives with the stream tests; tests are built from the repository only, never from an sdist.
    target_include_directories(obdeect_sources_test_${test_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../stream/cpp/tests)
    obdeect_sources_enable_warnings(obdeect_sources_test_${test_name})
    obdeect_add_test(${test_name} sources SOURCE test_${test_name}.cpp
                     COMMAND obdeect_sources_test_${test_name})
  endforeach()
  obdeect_check_test_sources(${CMAKE_CURRENT_SOURCE_DIR}/cpp/tests)
endif()
