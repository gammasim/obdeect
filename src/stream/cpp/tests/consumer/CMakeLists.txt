# A minimal consumer of the installed gstrickd package. Built by python/tests/test_cmake_package.py
# with -Dgstrickd_DIR=<directory of gstrickd-config.cmake>.
cmake_minimum_required(VERSION 3.20)
project(gstrickd_consumer LANGUAGES C CXX)

find_package(gstrickd CONFIG REQUIRED)

add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE gstrickd::gstrickd_core)

add_executable(consumer_c main_c.c)
target_link_libraries(consumer_c PRIVATE gstrickd::gstrickd_c)
