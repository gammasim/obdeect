# Registration of the native (ctest) tests of the obdeect components.
#
# Every test gets its component as a ctest label, a timeout, and the label "slow" when it
# takes longer than about 5 s (run `ctest -LE slow` for the fast set). A test with a
# source file is recorded, and obdeect_check_test_sources() fails the configure step when a
# cpp/tests/test_*.cpp file is not registered, so that no test is skipped silently.

# obdeect_add_test(<name> <component> [SOURCE <file.cpp>] [SLOW] COMMAND <command> [<args>...])
# SOURCE is the test source (a file name in cpp/tests) that the test runs; it is recorded for
# obdeect_check_test_sources().
function(obdeect_add_test name component)
  cmake_parse_arguments(ARG "SLOW" "SOURCE" "COMMAND" ${ARGN})
  if(NOT ARG_COMMAND)
    message(FATAL_ERROR "obdeect_add_test(${name}): COMMAND is required")
  endif()
  add_test(NAME ${name} COMMAND ${ARG_COMMAND})
  # Timeouts in seconds: 300 for a test, 900 for a slow one (a sanitizer build is several times slower).
  if(ARG_SLOW)
    set(labels "${component};slow")
    set(timeout 900)
  else()
    set(labels "${component}")
    set(timeout 300)
  endif()
  set_tests_properties(${name} PROPERTIES LABELS "${labels}" TIMEOUT ${timeout})
  if(ARG_SOURCE)
    get_filename_component(source_name "${ARG_SOURCE}" NAME)
    set_property(DIRECTORY APPEND PROPERTY OBDEECT_REGISTERED_TEST_SOURCES "${source_name}")
  endif()
endfunction()

# Fails the configure step when a test_*.cpp file of the directory is not registered with
# obdeect_add_test(... SOURCE ...) in the calling directory.
function(obdeect_check_test_sources tests_dir)
  file(GLOB sources CONFIGURE_DEPENDS "${tests_dir}/test_*.cpp")
  get_property(registered DIRECTORY PROPERTY OBDEECT_REGISTERED_TEST_SOURCES)
  foreach(source IN LISTS sources)
    get_filename_component(source_name "${source}" NAME)
    if(NOT source_name IN_LIST registered)
      message(FATAL_ERROR
              "obdeect: ${source_name} is not registered as a test in ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt "
              "(add it to the test list there, with its source)")
    endif()
  endforeach()
endfunction()
