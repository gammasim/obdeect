cmake_minimum_required(VERSION 3.20)
# The version comes from the git tag (setuptools-scm), not from here.
project(obdeect_ray_tracing LANGUAGES CXX)

option(OBDEECT_RAY_TRACING_BUILD_TESTING "Build the C++ tests" ON)
if(DEFINED SKBUILD_PLATLIB_DIR)
  set(OBDEECT_RAY_TRACING_BUILD_PYTHON_API_DEFAULT ON)
else()
  set(OBDEECT_RAY_TRACING_BUILD_PYTHON_API_DEFAULT OFF)
endif()
option(OBDEECT_RAY_TRACING_BUILD_BENCH "Build the benchmark and replay tools (not run by ctest)" OFF)
option(OBDEECT_RAY_TRACING_BUILD_PYTHON_API "Build the Python extension" ${OBDEECT_RAY_TRACING_BUILD_PYTHON_API_DEFAULT})

# The stream component: the target of the combined build, else the stream sources next to
# this component (builds from the repository), else the copy bundled into the sdist
# (stream-src, staged at release time), else an installed package.
if(NOT TARGET obdeect::stream)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../stream/CMakeLists.txt)
    set(OBDEECT_STREAM_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../stream)
  elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/stream-src/CMakeLists.txt)
    set(OBDEECT_STREAM_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/stream-src)
  endif()
  if(OBDEECT_STREAM_SOURCE_DIR)
    set(OBDEECT_STREAM_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(OBDEECT_STREAM_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
    set(OBDEECT_STREAM_BUILD_CAPI OFF CACHE BOOL "" FORCE)
    set(OBDEECT_STREAM_BUILD_PYTHON_API OFF CACHE BOOL "" FORCE)
    set(OBDEECT_STREAM_INSTALL OFF CACHE BOOL "" FORCE)
    add_subdirectory(${OBDEECT_STREAM_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/stream)
    message(STATUS "obdeect ray_tracing: using the stream sources in ${OBDEECT_STREAM_SOURCE_DIR}")
  else()
    find_package(Python 3.14 COMPONENTS Interpreter)
    if(Python_Interpreter_FOUND AND NOT obdeect_stream_DIR)
      execute_process(COMMAND "${Python_EXECUTABLE}" -m obdeect.stream --cmake-dir
                      OUTPUT_VARIABLE obdeect_stream_DIR OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    endif()
    find_package(obdeect_stream CONFIG REQUIRED)
    message(STATUS "obdeect ray_tracing: using the installed stream package ${obdeect_stream_VERSION} (${obdeect_stream_DIR})")
  endif()
endif()

# Header-only, standard-library-only core.
add_library(obdeect_ray_tracing_core INTERFACE)
target_include_directories(obdeect_ray_tracing_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include)
target_compile_features(obdeect_ray_tracing_core INTERFACE cxx_std_20)
# Philox core shared with the other components.
target_link_libraries(obdeect_ray_tracing_core INTERFACE obdeect::stream)
# Identical floating-point results for any thread count and block size on one
# platform: no implicit fused multiply-add contraction. Transcendental functions
# come from the platform's math library, so results may differ between
# platforms.
# -fno-math-errno lets GCC inline sqrt without the errno branch and libm call; it
# changes no result (the tracer never reads errno).
if(NOT MSVC)
  target_compile_options(obdeect_ray_tracing_core INTERFACE -ffp-contract=off -fno-math-errno)
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

function(obdeect_ray_tracing_enable_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /WX)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()
endfunction()

if(OBDEECT_RAY_TRACING_BUILD_PYTHON_API)
  find_package(Python 3.14 COMPONENTS Interpreter Development.Module REQUIRED)
  execute_process(COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
                  OUTPUT_VARIABLE nanobind_DIR OUTPUT_STRIP_TRAILING_WHITESPACE
                  COMMAND_ERROR_IS_FATAL ANY)
  find_package(nanobind CONFIG REQUIRED)
  nanobind_add_module(_core NB_STATIC cpp/bindings/core.cpp)
  target_link_libraries(_core PRIVATE obdeect_ray_tracing_core)
  obdeect_ray_tracing_enable_warnings(_core)
  install(TARGETS _core LIBRARY DESTINATION obdeect/ray_tracing)
endif()

if(OBDEECT_RAY_TRACING_BUILD_TESTING)
  enable_testing()
  foreach(name kernels trace model_json bounds)
    add_executable(test_${name} cpp/tests/test_${name}.cpp)
    target_link_libraries(test_${name} PRIVATE obdeect_ray_tracing_core)
    # check.hpp lives with the stream tests; tests are built from the repository only, never from an sdist.
    target_include_directories(test_${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../stream/cpp/tests)
    obdeect_ray_tracing_enable_warnings(test_${name})
    add_test(NAME ${name} COMMAND test_${name})
  endforeach()
  # Newton-step counts need the benchmark counters compiled in.
  add_executable(test_asphere_steps cpp/tests/test_asphere_steps.cpp)
  target_link_libraries(test_asphere_steps PRIVATE obdeect_ray_tracing_core)
  target_include_directories(test_asphere_steps PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../stream/cpp/tests)
  target_compile_definitions(test_asphere_steps PRIVATE OBDEECT_RAY_TRACING_BENCH_COUNTERS)
  obdeect_ray_tracing_enable_warnings(test_asphere_steps)
  add_test(NAME asphere_steps COMMAND test_asphere_steps)
  # The replay tool's limit comparison (header in cpp/bench, standard library only).
  add_executable(test_replay_compare cpp/tests/test_replay_compare.cpp)
  target_link_libraries(test_replay_compare PRIVATE obdeect_ray_tracing_core)
  target_include_directories(test_replay_compare PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../stream/cpp/tests)
  target_include_directories(test_replay_compare PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/bench)
  obdeect_ray_tracing_enable_warnings(test_replay_compare)
  add_test(NAME replay_compare COMMAND test_replay_compare)
endif()

if(OBDEECT_RAY_TRACING_BUILD_BENCH)
  foreach(name bench_trace replay)
    add_executable(${name} cpp/bench/${name}.cpp)
    target_link_libraries(${name} PRIVATE obdeect_ray_tracing_core)
    obdeect_ray_tracing_enable_warnings(${name})
  endforeach()
endif()
