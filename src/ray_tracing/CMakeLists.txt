cmake_minimum_required(VERSION 3.20)
# The version comes from the git tag (setuptools-scm), not from here.
project(obdeect_ray_tracing LANGUAGES CXX)

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

option(OBDEECT_RAY_TRACING_BUILD_TESTING "Build the C++ tests" ON)
option(OBDEECT_RAY_TRACING_BUILD_TOOLS "Build the native trace-stream program" ON)
if(DEFINED SKBUILD_PLATLIB_DIR)
  set(OBDEECT_RAY_TRACING_BUILD_PYTHON_API_DEFAULT ON)
else()
  set(OBDEECT_RAY_TRACING_BUILD_PYTHON_API_DEFAULT OFF)
endif()
option(OBDEECT_RAY_TRACING_BUILD_BENCH "Build the benchmark and replay tools (not run by ctest)" OFF)
option(OBDEECT_RAY_TRACING_BUILD_PYTHON_API "Build the Python extension" ${OBDEECT_RAY_TRACING_BUILD_PYTHON_API_DEFAULT})

# The stream component (see obdeect_locate_stream in cmake/ObdeectCommon.cmake).
obdeect_locate_stream(ray_tracing)

# Header-only, standard-library-only core.
find_package(Threads REQUIRED)
add_library(obdeect_ray_tracing_core INTERFACE)
target_include_directories(obdeect_ray_tracing_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include)
# Philox core shared with the other components.
target_link_libraries(obdeect_ray_tracing_core INTERFACE obdeect::stream Threads::Threads)
obdeect_set_core_options(obdeect_ray_tracing_core)
# -fno-math-errno lets GCC inline sqrt without the errno branch and libm call; it
# changes no result (the tracer never reads errno).
if(NOT MSVC)
  target_compile_options(obdeect_ray_tracing_core INTERFACE -fno-math-errno)
endif()

# Changes results at the rounding level. Use only on hardware with fused
# multiply-add; without it std::fma is emulated and very slow. Results differ
# from builds without it. Default builds, wheels and CI stay without FMA.
option(OBDEECT_RAY_TRACING_USE_FMA "Use std::fma in the aspheric polynomial evaluation" OFF)
if(OBDEECT_RAY_TRACING_USE_FMA)
  target_compile_definitions(obdeect_ray_tracing_core INTERFACE OBDEECT_RAY_TRACING_USE_FMA=1)
  # Without -mfma, x86-64 compilers call the math library for std::fma instead of inlining it.
  if(NOT MSVC AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|amd64")
    target_compile_options(obdeect_ray_tracing_core INTERFACE -mfma)
  endif()
endif()

if(OBDEECT_RAY_TRACING_BUILD_PYTHON_API)
  obdeect_add_python_module(ray_tracing obdeect_ray_tracing_core)
endif()

if(OBDEECT_RAY_TRACING_BUILD_TOOLS)
  add_executable(obdeect-ray-tracing-trace-stream cpp/tools/trace_stream.cpp)
  target_link_libraries(obdeect-ray-tracing-trace-stream PRIVATE obdeect_ray_tracing_core)
  obdeect_enable_warnings(obdeect-ray-tracing-trace-stream)
  if(OBDEECT_RAY_TRACING_BUILD_TESTING)
    obdeect_add_test(trace_stream_version ray_tracing COMMAND obdeect-ray-tracing-trace-stream --version)
    set_tests_properties(trace_stream_version PROPERTIES PASS_REGULAR_EXPRESSION
      "^obdeect ray-tracing trace-stream [0-9]+\\.[0-9]+\\.[0-9]+ \\(photon-stream format [0-9]+\\.[0-9]+\\)")
  endif()
  # Wheels: inside the package, run through `obdeect ray-tracing trace-stream`.
  if(DEFINED SKBUILD_PLATLIB_DIR)
    install(TARGETS obdeect-ray-tracing-trace-stream RUNTIME DESTINATION obdeect/ray_tracing/bin)
  endif()
endif()

if(OBDEECT_RAY_TRACING_BUILD_TESTING)
  enable_testing()
  # Tests that take more than about 5 s (labelled slow, see cmake/ObdeectTests.cmake).
  set(OBDEECT_RAY_TRACING_SLOW_TESTS kernels bounds)
  foreach(name kernels trace model_json bounds model_hash sphere_rows window)
    add_executable(test_${name} cpp/tests/test_${name}.cpp)
    target_link_libraries(test_${name} PRIVATE obdeect_ray_tracing_core)
    # check.hpp lives with the stream tests; tests are built from the repository only, never from an sdist.
    target_include_directories(test_${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../stream/cpp/tests)
    obdeect_enable_warnings(test_${name})
    set(slow_arg)
    if(name IN_LIST OBDEECT_RAY_TRACING_SLOW_TESTS)
      set(slow_arg SLOW)
    endif()
    obdeect_add_test(${name} ray_tracing SOURCE test_${name}.cpp ${slow_arg} COMMAND test_${name})
  endforeach()
  # The synthetic optical-model fixtures of the Python tests (labelled, no real telescope data).
  target_compile_definitions(test_trace PRIVATE
    OBDEECT_RAY_TRACING_TEST_MODELS="${CMAKE_CURRENT_SOURCE_DIR}/python/tests/data/optical-models")
  # The replay tool's limit comparison (header in cpp/bench, standard library only).
  add_executable(test_replay_compare cpp/tests/test_replay_compare.cpp)
  target_link_libraries(test_replay_compare PRIVATE obdeect_ray_tracing_core)
  target_include_directories(test_replay_compare PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../stream/cpp/tests)
  target_include_directories(test_replay_compare PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/bench)
  obdeect_enable_warnings(test_replay_compare)
  obdeect_add_test(replay_compare ray_tracing SOURCE test_replay_compare.cpp COMMAND test_replay_compare)
  obdeect_check_test_sources(${CMAKE_CURRENT_SOURCE_DIR}/cpp/tests)
endif()

if(OBDEECT_RAY_TRACING_BUILD_BENCH)
  foreach(name bench_trace replay)
    add_executable(${name} cpp/bench/${name}.cpp)
    target_link_libraries(${name} PRIVATE obdeect_ray_tracing_core)
    obdeect_enable_warnings(${name})
  endforeach()
endif()
