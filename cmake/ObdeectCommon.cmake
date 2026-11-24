# Build helpers shared by the obdeect components (src/<component>/CMakeLists.txt).
#
# A component includes this file through OBDEECT_SHARED_CMAKE_DIR: the repository's cmake/
# directory when it builds from the repository, or the copy staged as shared-cmake/ next to
# the component when it builds from an sdist (see .github/workflows/pypi.yml).
include_guard(GLOBAL)

include("${CMAKE_CURRENT_LIST_DIR}/ObdeectTests.cmake")

# Warnings are errors for the targets of a component (tests, tools, bindings).
function(obdeect_enable_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /WX)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()
endfunction()

# Options of a header-only core target: C++20 and no implicit fused multiply-add
# contraction, so results do not depend on the thread count or the block or chunk size on
# one platform. Transcendental functions come from the platform's math library, so results
# may differ between platforms.
function(obdeect_set_core_options target)
  target_compile_features(${target} INTERFACE cxx_std_20)
  if(NOT MSVC)
    target_compile_options(${target} INTERFACE -ffp-contract=off)
  endif()
endfunction()

# Makes obdeect::stream available to a component. Uses, in order: the stream target of the
# combined build; the stream sources next to the component (repository builds); the stream
# sources staged in the component sdist (stream-src/); an installed obdeect_stream package.
# Stream's own CMakeLists.txt does not use it.
macro(obdeect_locate_stream component)
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
      message(STATUS "obdeect ${component}: using the stream sources in ${OBDEECT_STREAM_SOURCE_DIR}")
    else()
      find_package(Python 3.14 COMPONENTS Interpreter)
      if(Python_Interpreter_FOUND AND NOT obdeect_stream_DIR)
        execute_process(COMMAND "${Python_EXECUTABLE}" -m obdeect.stream --cmake-dir
                        OUTPUT_VARIABLE obdeect_stream_DIR OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
      endif()
      find_package(obdeect_stream CONFIG REQUIRED)
      message(STATUS "obdeect ${component}: using the installed stream package ${obdeect_stream_VERSION} (${obdeect_stream_DIR})")
    endif()
  endif()
endmacro()

# The Python extension obdeect/<component>/_core, built from cpp/bindings/core.cpp of the
# calling directory and linked to the core target. nanobind is the only Python dependency of
# the C++ build and it is needed only here.
function(obdeect_add_python_module component core_target)
  find_package(Python 3.14 COMPONENTS Interpreter Development.Module REQUIRED)
  execute_process(COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
                  OUTPUT_VARIABLE nanobind_DIR OUTPUT_STRIP_TRAILING_WHITESPACE
                  COMMAND_ERROR_IS_FATAL ANY)
  find_package(nanobind CONFIG REQUIRED)
  nanobind_add_module(_core NB_STATIC cpp/bindings/core.cpp)
  target_link_libraries(_core PRIVATE ${core_target})
  obdeect_enable_warnings(_core)
  install(TARGETS _core LIBRARY DESTINATION obdeect/${component})
endfunction()
